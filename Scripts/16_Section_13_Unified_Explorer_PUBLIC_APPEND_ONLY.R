# ==============================================================================
# SCRIPT 16 — APPEND-ONLY SECTION 13: UNIFIED INTERACTIVE UMAP EXPLORER
# ==============================================================================
# PURPOSE
#   Provide code for a GitHub Pages-compatible explorer to append to the END of
#   the ORIGINAL Script 16, without changing existing NMF H processing, UMAP
#   fitting, statistics, CSVs, PDFs, PNGs, or individual Plotly HTML files.
#
# WHERE TO PASTE
#   Paste THIS ENTIRE SECTION AFTER the existing:
#     # END OF SCRIPT 16
#   The earlier script will run fully; only then does this section execute.
#
# WHAT THE DASHBOARD DOES
#   * Dropdown to color all patients by any approved covariate.
#   * Dropdown to toggle between raw H and normalized H UMAP coordinates.
#   * Dropdown to highlight any NMF group while dimming other patients.
#   * Hover tooltips, zoom/pan, legend toggles, and Plotly PNG export.
#   * Never recomputes NMF or UMAP. Uses SAVED primary embeddings.
#
# METHODS / INTERPRETATION
#   Recoloring or highlighting does not change the UMAP coordinates. Raw and
#   normalized H are independently fitted embeddings; their axes and distances
#   should not be compared as if on a shared coordinate system. The explorer
#   is for descriptive visualization, not statistical hypothesis testing.
#
# ==============================================================================

# ==============================================================================
# SECTION 13.01 — EXPLORER-ONLY SETTINGS (DOES NOT MODIFY CONFIG)
# ==============================================================================
EXPLORER_CONFIG <- list(
  ENABLED = FALSE,  # Keep FALSE until public dataset/derived-data release is approved.
  PUBLIC_RELEASE_APPROVED = FALSE,
  # EXPLICIT allowlist only; do not export all CONFIG$PLOT_COVARIATES.
  PUBLIC_COLOR_COVARIATES = character(0),
  # These extra variables appear on hover and are embedded in the HTML.
  # Only add columns approved for public redistribution.
  PUBLIC_TOOLTIP_COVARIATES = character(0),
  # TRUE = single portable HTML for GitHub Pages; requires Pandoc.
  SELF_CONTAINED = TRUE
)

# ==============================================================================
# SECTION 13.02 — LOAD SAVED RESULTS; SAFELY APPEND TO EXISTING LOG
# ==============================================================================
# The main Script 16 closes its log connection before this new section starts.
# We therefore append directly to its existing text log (no reused connection).
EXPLORER_DIR <- if (exists("OUTPUT_DIR", inherits = TRUE)) {
  OUTPUT_DIR
} else {
  file.path(CONFIG$ROOT_DIR, CONFIG$OUTPUT_DIR)
}
if (!dir.exists(EXPLORER_DIR)) {
  stop("Section 13 requires the completed Script 16 output directory: ",
       EXPLORER_DIR)
}
EXPLORER_LOG <- file.path(
  EXPLORER_DIR, "16_NMF_H_UMAP_Interactive_Exploration_Log.txt"
)
explorer_log <- function(..., level = "INFO") {
  message_text <- paste(..., collapse = "")
  line <- sprintf(
    "[%s] [%s] [SECTION 13 - EXPLORER] %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, message_text
  )
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = EXPLORER_LOG, append = TRUE)
  invisible(line)
}
explorer_record_output <- function(path) {
  explorer_log("Saved: ", normalizePath(path, mustWork = TRUE))
}

# This saved object contains raw/normalized annotated UMAP coordinates from
# the ORIGINAL analysis and does not require new UMAP computation.
EXPLORER_SOURCE_RDS <- file.path(
  EXPLORER_DIR, "Script_16_Complete_Analysis_Objects.rds"
)

# ==============================================================================
# SECTION 13.03 — STRICT COVARIATE TYPING AND COLOR PALETTES
# ==============================================================================
# These helpers intentionally live inside the append-only section: they do not
# rely on the original log connection and do not overwrite old helper functions.
explorer_make_field <- function(df, field) {
  x <- df[[field]]
  if (field %in% CONFIG$BINARY_LESION_COLS) {
    raw <- tolower(trimws(as.character(x)))
    present <- c("1", "true", "t", "yes", "y", "positive", "pos", "+", "present")
    absent <- c("0", "false", "f", "no", "n", "negative", "neg", "-", "absent")
    out <- rep("Unknown", length(raw))
    out[!is.na(raw) & raw %in% present] <- "Present"
    out[!is.na(raw) & raw %in% absent] <- "Absent"
    unexpected <- unique(raw[!is.na(raw) & nzchar(raw) &
                               !raw %in% c(present, absent, "na", "n/a",
                                           "unknown", "missing", "not tested",
                                           "not available")])
    if (length(unexpected)) {
      explorer_log("Binary field '", field, "' had unrecognized values: ",
                   paste(head(unexpected, 10L), collapse = ", "),
                   "; encoded as Unknown", level = "WARNING")
    }
    return(factor(out, levels = c("Absent", "Present", "Unknown")))
  }
  if (field %in% CONFIG$FORCE_CATEGORICAL) {
    out <- as.character(x)
    out[is.na(out) | !nzchar(trimws(out))] <- "Unknown"
    return(factor(out))
  }
  if (field %in% CONFIG$FORCE_CONTINUOUS || is.numeric(x)) {
    nums <- suppressWarnings(as.numeric(as.character(x)))
    if (field %in% CONFIG$FORCE_CONTINUOUS &&
        any(!is.na(x) & is.na(nums))) {
      stop("Explorer forced-continuous column contains nonnumeric values: ", field)
    }
    return(nums)
  }
  out <- as.character(x)
  out[is.na(out) | !nzchar(trimws(out))] <- "Unknown"
  factor(out)
}
explorer_palette <- function(field, values) {
  if (field %in% CONFIG$BINARY_LESION_COLS) {
    return(CONFIG$LESION_COLORS)
  }
  chosen <- CONFIG$CATEGORICAL_PALETTES[[field]]
  if (is.null(chosen)) chosen <- character(0)
  categories <- sort(unique(as.character(values)))
  missing <- setdiff(categories, names(chosen))
  if (length(missing)) {
    additional <- grDevices::hcl(
      h = seq(15, 375, length.out = length(missing) + 1L)[
        seq_len(length(missing))
      ], c = 70, l = 55
    )
    names(additional) <- missing
    chosen <- c(chosen, additional)
  }
  chosen
}

# ==============================================================================
# SECTION 13.04 — DASHBOARD CONSTRUCTION
# ==============================================================================
# Relies only on the saved Script 16 embeddings and the functions above.
# The browser-side JavaScript is embedded through htmlwidgets::onRender().
# GitHub Pages is sufficient: there is no Shiny server or other R backend.
make_unified_explorer <- function(embeddings) {
  if (!isTRUE(EXPLORER_CONFIG$ENABLED)) {
    explorer_log("Unified explorer disabled by configuration")
    return(invisible(NULL))
  }
  if (!all(c("raw", "normalized") %in% names(embeddings))) {
    stop("Unified explorer requires both raw and normalized H embeddings. ",
         "Set REPRESENTATIONS = c('raw', 'normalized').")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Unified explorer requires package 'jsonlite'. Install it first.")
  }

  raw <- embeddings[["raw"]]
  normalized <- embeddings[["normalized"]]
  id <- CONFIG$PATIENT_ID_COL
  if (anyDuplicated(raw[[id]]) || anyDuplicated(normalized[[id]])) {
    stop("Explorer has duplicate patient identifiers")
  }
  j <- match(raw[[id]], normalized[[id]])
  if (anyNA(j) || nrow(raw) != nrow(normalized)) {
    stop("Raw and normalized embeddings contain different patients")
  }
  normalized <- normalized[j, , drop = FALSE]
  if (!identical(as.character(raw[[id]]), as.character(normalized[[id]]))) {
    stop("Failed to align H representations by patient identifier")
  }

  approved_colors <- EXPLORER_CONFIG$PUBLIC_COLOR_COVARIATES
  if (is.null(approved_colors)) approved_colors <- character(0)
  approved_colors <- unique(as.character(approved_colors))
  available_colors <- intersect(approved_colors, names(raw))
  missing_colors <- setdiff(approved_colors, names(raw))
  if (length(missing_colors)) {
    explorer_log("Explorer requested nonexistent color columns: ",
             paste(missing_colors, collapse = ", "), level = "WARNING")
  }
  if (!length(available_colors)) stop("No approved explorer coloring fields found")
  if (!"NMF_Group" %in% names(raw)) stop("NMF_Group required for group highlighting")
  if (!all(available_colors %in% names(normalized))) stop("Covariate mismatch between representations")

  approved_hover <- unique(as.character(EXPLORER_CONFIG$PUBLIC_TOOLTIP_COVARIATES))
  approved_hover <- intersect(approved_hover, names(raw))
  # Only needed columns are serialized into the public-facing HTML document.
  data_cols <- unique(c("NMF_Group", available_colors, approved_hover))
  # Direct patient identifiers are not permitted in the public export.
  forbidden_fields <- c(id, "USI", "Public_ID", "OS.time..days.", "OS.event.ID")
  if (any(c(approved_colors, approved_hover) %in% forbidden_fields)) {
    stop("Public explorer requested a direct identifier or disallowed outcome field")
  }
  export <- data.frame(
    Public_ID = sprintf("AML_%04d", seq_len(nrow(raw))),
    RawX = as.numeric(raw$UMAP1), RawY = as.numeric(raw$UMAP2),
    NormX = as.numeric(normalized$UMAP1),
    NormY = as.numeric(normalized$UMAP2), stringsAsFactors = FALSE
  )
  for (field in data_cols) {
    export[[field]] <- explorer_make_field(raw, field)
    if (is.factor(export[[field]])) export[[field]] <- as.character(export[[field]])
  }
  if (any(!is.finite(unlist(export[c("RawX", "RawY", "NormX", "NormY")])))) {
    stop("Explorer UMAP coordinates contain non-finite values")
  }

  # Determine categorical versus numeric presentation and its stable palette.
  descriptors <- lapply(available_colors, function(field) {
    x <- export[[field]]
    numeric <- is.numeric(x)
    list(name = field,
         numeric = numeric,
         palette = if (numeric) NULL else as.list(explorer_palette(field, x)))
  })
  named_palette <- list()
  for (d in descriptors) if (!d$numeric) named_palette[[d$name]] <- d$palette
  # Explicit record encoding means the JavaScript receives an array of objects.
  records <- lapply(seq_len(nrow(export)), function(i) as.list(export[i, , drop = FALSE]))
  payload <- list(
    rows = records,
    colorFields = available_colors,
    hoverFields = approved_hover,
    numericFields = vapply(descriptors, function(d) d$numeric, logical(1)),
    palettes = named_palette,
    groups = sort(unique(as.character(export$NMF_Group))),
    defaultColor = available_colors[[1L]],
    low = CONFIG$CONTINUOUS_LOW_COLOR,
    mid = CONFIG$CONTINUOUS_MID_COLOR,
    high = CONFIG$CONTINUOUS_HIGH_COLOR,
    title = CONFIG$PLOT_TITLE_PREFIX
  )

  # Minimal invisible seed trace loads the bundled plotly.js runtime. After the
  # widget renders, the controls update the figure entirely in the browser.
  widget <- plotly::plot_ly(x = c(0), y = c(0), type = "scatter",
                            mode = "markers", showlegend = FALSE,
                            marker = list(opacity = 0), hoverinfo = "skip")
  widget <- plotly::layout(widget,
                           title = list(text = "Interactive NMF H UMAP Explorer"),
                           xaxis = list(title = "UMAP 1", zeroline = FALSE),
                           yaxis = list(title = "UMAP 2", zeroline = FALSE),
                           margin = list(t = 55, l = 65, r = 30, b = 65))

  # Every DOM selector and every point is created from a controlled local
  # dataset. No external scripts, Shiny services or per-covariate HTML files.
  javascript <- '
function(el, x, info) {
  if (el.dataset.nmfExplorerInit === "1") return;
  el.dataset.nmfExplorerInit = "1";
  const rows = info.rows || [];
  const holder = document.createElement("div");
  holder.style.cssText = "display:flex;flex-wrap:wrap;gap:12px;align-items:end;padding:14px;background:#f5f7fb;border:1px solid #dfe5ee;border-radius:9px;font:14px Arial,sans-serif;color:#182236;margin-bottom:12px";
  const createSelect = function(title, values, display) {
    const lab = document.createElement("label");
    lab.style.cssText = "display:flex;flex-direction:column;gap:5px;font-weight:600;min-width:170px";
    const caption = document.createElement("span"); caption.textContent = title;
    const sel = document.createElement("select");
    sel.style.cssText = "border:1px solid #b7c3d5;border-radius:6px;padding:8px;background:white;color:#17223a;width:100%;max-width:300px";
    values.forEach(function(value) {
      const option = document.createElement("option");
      option.value = String(value);
      option.textContent = display ? display(value) : String(value);
      sel.appendChild(option);
    });
    lab.appendChild(caption); lab.appendChild(sel); holder.appendChild(lab);
    return sel;
  };
  const colorSelect = createSelect("Color by", info.colorFields);
  colorSelect.value = info.defaultColor;
  const representationSelect = createSelect("H representation", ["normalized", "raw"],
    v => v === "normalized" ? "Normalized H" : "Raw H");
  const groupSelect = createSelect("Highlight NMF group", ["all"].concat(info.groups),
    v => v === "all" ? "All patients" : v);
  const count = document.createElement("span");
  count.style.cssText = "font-weight:500;font-size:12px;color:#4e596b;padding-bottom:8px";
  holder.appendChild(count);
  const scaleBar = document.createElement("div");
  scaleBar.style.cssText = "display:none;align-self:end;min-width:170px;font-size:11px;color:#344054";
  holder.appendChild(scaleBar);
  el.parentNode.insertBefore(holder, el);
  function escapeHTML(value) {
    return String(value == null ? "Unknown" : value)
      .replace(/&/g,"&amp;").replace(/</g,"&lt;")
      .replace(/>/g,"&gt;").replace(/\"/g,"&quot;").replace(/\x27/g,"&#39;");
  }
  function hexRGB(hex) {
    const h = String(hex).replace("#", "");
    return [parseInt(h.substring(0,2),16),parseInt(h.substring(2,4),16),parseInt(h.substring(4,6),16)];
  }
  function blend(a,b,t) {
    const c=hexRGB(a),d=hexRGB(b),q=Math.max(0,Math.min(1,t));
    return "rgb("+c.map((v,i)=>Math.round(v+(d[i]-v)*q)).join(",")+")";
  }
  function numericColor(v, min, max) {
    if (v===null || v===undefined || v==="" || !Number.isFinite(Number(v))) return "#BDBDBD";
    const frac = max===min ? 0.5 : (Number(v)-min)/(max-min);
    return frac<0.5 ? blend(info.low,info.mid,frac*2) : blend(info.mid,info.high,(frac-.5)*2);
  }
  function hover(row,variable) {
    const fields = Array.from(new Set(["Public_ID", variable].concat(info.hoverFields)));
    return fields.map(function(field) {
      const value = row[field] === null || row[field] === undefined || row[field] === "" ? "Unknown" : row[field];
      return "<b>"+escapeHTML(field)+":</b> "+escapeHTML(value);
    }).join("<br>");
  }
  function makeTrace(subset, label, color, rep, variable, numeric, opacity) {
    const xKey = rep === "raw" ? "RawX" : "NormX";
    const yKey = rep === "raw" ? "RawY" : "NormY";
    return {type:"scattergl",mode:"markers",name:label,
      x:subset.map(r=>r[xKey]), y:subset.map(r=>r[yKey]),
      text:subset.map(r=>hover(r,variable)),
      hovertemplate:"%{text}<extra></extra>",
      marker:{color:color,size:7,opacity:opacity||0.87},
      showlegend:!!label};
  }
  function draw() {
    const rep=representationSelect.value,variable=colorSelect.value,selected=groupSelect.value;
    const numeric = !!info.numericFields[info.colorFields.indexOf(variable)];
    const focus = selected === "all" ? rows : rows.filter(r=>String(r.NMF_Group)===selected);
    const others = selected === "all" ? [] : rows.filter(r=>String(r.NMF_Group)!==selected);
    const traces=[];
    if (others.length) traces.push(makeTrace(others,"Other patients","#cccccc",rep,variable,false,0.15));
    scaleBar.style.display = numeric ? "block" : "none";
    if (numeric) {
      const valid=rows.map(r=>r[variable] == null ? NaN : Number(r[variable])).filter(Number.isFinite);
      const min=valid.length?Math.min(...valid):0,max=valid.length?Math.max(...valid):1;
      scaleBar.innerHTML = "";
      const gradient=document.createElement("div");
      gradient.style.cssText="height:10px;border-radius:5px;background:linear-gradient(90deg,"+info.low+","+info.mid+","+info.high+")";
      const limits=document.createElement("div");
      limits.style.cssText="display:flex;justify-content:space-between;margin-top:4px";
      const minLabel=document.createElement("span"); minLabel.textContent=Number(min).toPrecision(3);
      const maxLabel=document.createElement("span"); maxLabel.textContent=Number(max).toPrecision(3);
      limits.appendChild(minLabel); limits.appendChild(maxLabel);
      scaleBar.appendChild(gradient); scaleBar.appendChild(limits);
      const colors=focus.map(r=>numericColor(r[variable],min,max));
      traces.push(makeTrace(focus,null,colors,rep,variable,true,0.9));
      traces[traces.length-1].marker.colorscale=[[0,info.low],[0.5,info.mid],[1,info.high]];
      traces[traces.length-1].marker.showscale=false;
    } else {
      const categories=Array.from(new Set(rows.map(r=>String(r[variable]==null?"Unknown":r[variable])))).sort();
      const palette=info.palettes[variable]||{};
      categories.forEach(function(category) {
        const selectedRows=focus.filter(r=>String(r[variable]==null?"Unknown":r[variable])===category);
        if (!selectedRows.length) return;
        traces.push(makeTrace(selectedRows,category,palette[category]||"#888888",rep,variable,false,0.87));
      });
    }
    count.textContent = focus.length + " of " + rows.length + " patients highlighted";
    const layout={title:{text:info.title+" — "+(rep==="raw"?"Raw H":"Normalized H")},
      xaxis:{title:"UMAP 1",zeroline:false},yaxis:{title:"UMAP 2",zeroline:false,scaleanchor:"x",scaleratio:1},
      hovermode:"closest",legend:{orientation:"v",x:1.01,y:1},
      margin:{l:65,r:50,b:70,t:65},paper_bgcolor:"white",plot_bgcolor:"white"};
    Plotly.react(el,traces,layout,{responsive:true,displaylogo:false,
      toImageButtonOptions:{format:"png",filename:"NMF_UMAP_"+rep+"_"+variable,scale:3}});
  }
  colorSelect.addEventListener("change",draw);
  representationSelect.addEventListener("change",draw);
  groupSelect.addEventListener("change",draw);
  draw();
}
'
  widget <- htmlwidgets::onRender(widget, javascript, data = payload)
  destination <- file.path(EXPLORER_DIR, "PUBLIC_RELEASE_REVIEW_REQUIRED",
                           "Unified_UMAP_Explorer_PUBLIC_CANDIDATE.html")
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  # A truly stand-alone Github Pages page should not require a separate
  # dependencies directory. Raise an error rather than silently publishing a
  # broken page if Pandoc/selfcontained conversion is unavailable.
  htmlwidgets::saveWidget(widget, destination,
                          selfcontained = isTRUE(EXPLORER_CONFIG$SELF_CONTAINED),
                          title = "Pediatric AML NMF — Interactive UMAP Explorer")
  explorer_record_output(destination)
  explorer_log("Unified explorer created: ", destination)
  explorer_log("Embedded anonymized rows: ", nrow(export))
  explorer_log("Public color variables: ", paste(available_colors, collapse = ", "))
  explorer_log("Public tooltip-only variables: ", paste(setdiff(approved_hover, available_colors), collapse = ", "))
  explorer_log("IMPORTANT: Review the HTML and embedded patient-level attributes ",
           "before public GitHub publication.", level = "WARNING")
  invisible(destination)
}

# ==============================================================================
# SECTION 13.05 — EXECUTE EXPLORER ONLY; RECORD SUCCESS OR ERROR
# ==============================================================================
if (isTRUE(EXPLORER_CONFIG$ENABLED)) {
  if (!isTRUE(EXPLORER_CONFIG$PUBLIC_RELEASE_APPROVED)) {
    stop("Public explorer export requires explicit release approval. ",
         "Keep ENABLED=FALSE for publicly distributed code.")
  }
  if (!length(EXPLORER_CONFIG$PUBLIC_COLOR_COVARIATES)) {
    stop("Set an explicit nonempty PUBLIC_COLOR_COVARIATES allowlist")
  }
  explorer_log("Starting public-candidate explorer (review required before release)")
  explorer_failure <- tryCatch({
    for (pkg in c("plotly", "htmlwidgets", "jsonlite")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Missing package '", pkg, "'. Install with install.packages('", pkg, "').")
      }
    }
    if (!file.exists(EXPLORER_SOURCE_RDS)) {
      stop("Original Script 16 analysis object not found: ", EXPLORER_SOURCE_RDS)
    }
    source_obj <- readRDS(EXPLORER_SOURCE_RDS)
    if (is.null(source_obj$embeddings) ||
        !all(c("raw", "normalized") %in% names(source_obj$embeddings))) {
      stop("Saved Script 16 object lacks raw and/or normalized UMAP embeddings")
    }
    destination <- make_unified_explorer(source_obj$embeddings)
    # Add just the NEW deliverable to the existing manifest; do not overwrite it.
    manifest <- file.path(EXPLORER_DIR, "Output_Manifest.txt")
    cat(normalizePath(destination, mustWork = TRUE), "\n",
        sep = "", file = manifest, append = TRUE)
    explorer_log("Public CANDIDATE completed. No permission to publish is implied.")
    NULL
  }, error = function(e) {
    explorer_log("FATAL ERROR: ", conditionMessage(e), level = "ERROR")
    e
  })
  if (inherits(explorer_failure, "error")) stop(explorer_failure)
}
# ==============================================================================
# END OF APPEND-ONLY SECTION 13 — ORIGINAL SCRIPT 16 REMAINS UNCHANGED
# ==============================================================================
