# ==============================================================================
# SCRIPT 16 — SECTION 14: PUBLIC GITHUB UMAP EXPLORER EXPORTER (STANDALONE)
# ==============================================================================
# SCIENTIFIC PURPOSE
#   Build one interactive GitHub Pages-compatible Plotly HTML directly from
#   the RAW and NORMALIZED H-matrix UMAP coordinates saved by Script 16.
#   This script does NOT recalculate NMF, H normalization, or UMAP.
#
# PUBLIC-DATA DESIGN
#   - No source patient identifiers (USI) in the published records.
#   - Assign newly randomized, generic AML_0001 ... display IDs.
#   - NO survival, vital status, event, death, follow-up or time-to-event
#     variables may be selected for coloring or tooltips.
#   - Only explicitly allowlisted annotations are serialized into HTML.
#   - Original USI-to-display-label mapping is NOT written to disk.
#   - Private audit/log saved adjacent to the original analysis, not to GitHub.
#
# IMPORTANT PRIVACY / DATA-USE LIMITATION
#   Display ID replacement is PSEUDONYMIZATION, not complete anonymization.
#   Patient-specific UMAP positions and combinations of rare molecular or
#   clinical annotations remain recoverable from the HTML source.
#
# EXPECTED INPUT
#   Script_16_Complete_Analysis_Objects.rds (created by original Script 16)
#   with $embeddings$raw and $embeddings$normalized data frames containing
#   USI, UMAP1, UMAP2, NMF_Group, and covariate columns.
#
# EXPECTED OUTPUTS
#   Unified_UMAP_Explorer_PUBLIC.html       
#   Public_Explorer_Field_Audit_PRIVATE.csv 
#   Public_Explorer_Export_Log_PRIVATE.txt  
#
# DEPENDENCIES: plotly, htmlwidgets, jsonlite; Pandoc for self-contained HTML
# ============================================================================== 

# ==============================================================================
# SECTION 14.01 — USER CONFIGURATION: ONLY EDIT THIS AREA
# ==============================================================================
# Environment variable is preferred, but setting this path directly works too.
# Point to the existing Script 16 output directory; no model input files needed.
EXPLORER_CONFIG <- list(
  SCRIPT16_OUTPUT_DIR = Sys.getenv("SCRIPT16_OUTPUT_DIR", unset = ""),
  # Example alternative (uncomment and customize ONLY on your private laptop):
  # SCRIPT16_OUTPUT_DIR = "/path/to/16_NMF_H_UMAP_Interactive_Exploration",
  SOURCE_RDS = "Script_16_Complete_Analysis_Objects.rds",
  PATIENT_ID_COL = "USI",
  SELF_CONTAINED = TRUE,
  PUBLIC_ID_SEED = 190L,
  # Your complete intended PUBLIC covariate menu; unavailable ones are logged.
  # Keep this an explicit allowlist rather than exporting all clinical fields.
  PUBLIC_COLOR_COVARIATES = c(
    "NMF_Group", "Sex", "Risk.group.classification..by.protocol.1031",
    "Margin", "Margin_Relative", "Age.Category", "AgeGroup",
    "M6_AML", "M7_AML", "FLT3.ITD.positive.", "WT1.mutation.",
    "NPM.mutation.", "CEBPA.mutation.", "RAS.Gene", "Mutations.Category",
    "Primary.CNV", "Primary.Fusion", "Cytogenetic.Category.1",
    "Cytogenetic.Category.2", "t.9.11..p22.q23.", "t.8.21",
    "del5q", "inv.16.", "trisomy.21", "MLL", "Minus.Y",
    "Minus.X", "Primary.Cytogenetic.Code"
  ),
  # Tooltips show ONLY these details plus the actively selected color variable
  # and generic Public_ID; NMF group is also needed for highlighting.
  PUBLIC_TOOLTIP_COVARIATES = c(
    "NMF_Group", "Sex", "Risk.group.classification..by.protocol.1031",
    "Margin", "Margin_Relative"
  ),
  BINARY_LESION_COLS = c(
    "M6_AML", "M7_AML", "FLT3.ITD.positive.", "WT1.mutation.",
    "NPM.mutation.", "CEBPA.mutation.", "t.9.11..p22.q23.",
    "t.8.21", "del5q", "inv.16.", "trisomy.21", "MLL",
    "Minus.Y", "Minus.X"
  ),
  FORCE_CONTINUOUS = c("Margin", "Margin_Relative"),
  FORCE_CATEGORICAL = c(
    "NMF_Group", "Sex", "Risk.group.classification..by.protocol.1031",
    "Age.Category", "AgeGroup", "RAS.Gene", "Mutations.Category",
    "Primary.CNV", "Primary.Fusion", "Cytogenetic.Category.1",
    "Cytogenetic.Category.2", "Primary.Cytogenetic.Code"
  ),
  CATEGORICAL_PALETTES = list(
    Sex = c(Male="#2B6CB0", Female="#DD6B98", Unknown="#999999"),
    Risk.group.classification..by.protocol.1031 = c(
      Low="#239A3B", High="#D33F49",
      Unknown="#999999")
  ),
  LESION_COLORS = c(Present="#D62728", Absent="#8C8C8C", Unknown="#D8D8D8"),
  CONTINUOUS_LOW_COLOR = "#2B6CB0",
  CONTINUOUS_MID_COLOR = "#F7F7F7",
  CONTINUOUS_HIGH_COLOR = "#C53030",
  PLOT_TITLE_PREFIX = "Pediatric AML NMF H UMAP",
  # Explicit manual acknowledgement. FALSE refuses to export real patient rows.
  REDISTRIBUTION_REVIEW_COMPLETED = FALSE
)

# Compatibility aliases used in the reusable explorer helper functions below.
CONFIG <- list(
  PATIENT_ID_COL = EXPLORER_CONFIG$PATIENT_ID_COL,
  BINARY_LESION_COLS = EXPLORER_CONFIG$BINARY_LESION_COLS,
  FORCE_CATEGORICAL = EXPLORER_CONFIG$FORCE_CATEGORICAL,
  FORCE_CONTINUOUS = EXPLORER_CONFIG$FORCE_CONTINUOUS,
  CATEGORICAL_PALETTES = EXPLORER_CONFIG$CATEGORICAL_PALETTES,
  LESION_COLORS = EXPLORER_CONFIG$LESION_COLORS,
  CONTINUOUS_LOW_COLOR = EXPLORER_CONFIG$CONTINUOUS_LOW_COLOR,
  CONTINUOUS_MID_COLOR = EXPLORER_CONFIG$CONTINUOUS_MID_COLOR,
  CONTINUOUS_HIGH_COLOR = EXPLORER_CONFIG$CONTINUOUS_HIGH_COLOR,
  PLOT_TITLE_PREFIX = EXPLORER_CONFIG$PLOT_TITLE_PREFIX
)

# ==============================================================================
# SECTION 14.02 — SAFETY: PATHS, ALLOWLIST, OUTCOME/IDENTIFIER EXCLUSION
# ==============================================================================
BLOCKED_FIELD_REGEX <- paste0(
  "(^|[._])(",
  "usi|patient.?id|sample.?id|subject.?id|case.?id|barcode|",
  "surv|survival|overall.?survival|progression|relapse|",
  "death|died|deceased|mortality|vital|event|follow.?up|",
  "os|efs|pfs|dfs|rfs|days.?to|time.?to|",
  "date.?of|birth|dob|contact|address|phone|email|",
  "accession|uuid|file.?id",
  ")([._]|$)"
)
assert_safe_fields <- function(fields, purpose) {
  fields <- as.character(fields)
  # Public_ID is a local generated label, NOT source identifying data.
  check <- fields[fields != "Public_ID"]
  forbidden <- check[grepl(BLOCKED_FIELD_REGEX, check, ignore.case=TRUE, perl=TRUE)]
  if (length(forbidden)) {
    stop("Unsafe ", purpose, " field(s) are prohibited from HTML: ",
         paste(unique(forbidden), collapse=", "))
  }
  invisible(TRUE)
}

if (!isTRUE(EXPLORER_CONFIG$REDISTRIBUTION_REVIEW_COMPLETED)) {
  stop("Export intentionally disabled: first verify that ALL source datasets, ",
       "derived UMAP coordinates, and requested individual-level annotations ",
       "may be redistributed publicly. After completing this review, set ",
       "REDISTRIBUTION_REVIEW_COMPLETED = TRUE. This is a user acknowledgement, ",
       "NOT an automated determination of legal/data-use permission.")
}
if (!nzchar(EXPLORER_CONFIG$SCRIPT16_OUTPUT_DIR)) {
  stop("Set SCRIPT16_OUTPUT_DIR environment variable or configuration to the ",
       "existing private Script 16 results folder.")
}
EXPLORER_DIR <- normalizePath(EXPLORER_CONFIG$SCRIPT16_OUTPUT_DIR, mustWork=TRUE)
EXPLORER_SOURCE_RDS <- file.path(EXPLORER_DIR, EXPLORER_CONFIG$SOURCE_RDS)
EXPLORER_LOG <- file.path(EXPLORER_DIR, "Public_Explorer_Export_Log_PRIVATE.txt")
explorer_log <- function(..., level="INFO") {
  ln <- sprintf("[%s] [%s] [PUBLIC HTML EXPORT] %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste(...,collapse=""))
  cat(ln, "\n", sep="")
  cat(ln, "\n", sep="", file=EXPLORER_LOG, append=TRUE)
  invisible(ln)
}
explorer_record_output <- function(path) {
  explorer_log("Saved: ", basename(path))
}
assert_safe_fields(EXPLORER_CONFIG$PUBLIC_COLOR_COVARIATES, "configured color")
assert_safe_fields(EXPLORER_CONFIG$PUBLIC_TOOLTIP_COVARIATES, "configured tooltip")

# ==============================================================================
# SECTION 14.03 — PACKAGES AND SAVED OBJECT VALIDATION
# ==============================================================================
for (pkg in c("plotly", "htmlwidgets", "jsonlite")) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    stop("Missing R package: ", pkg, ". Run install.packages('",pkg,"').")
  }
}
if (!file.exists(EXPLORER_SOURCE_RDS)) stop("Missing Script 16 object: ",EXPLORER_SOURCE_RDS)
source_obj <- readRDS(EXPLORER_SOURCE_RDS)
if (is.null(source_obj$embeddings) ||
    !all(c("raw", "normalized") %in% names(source_obj$embeddings))) {
  stop("Saved object does not contain embeddings$raw and embeddings$normalized")
}
for (nm in c("raw", "normalized")) {
  dat <- source_obj$embeddings[[nm]]
  required <- c(CONFIG$PATIENT_ID_COL, "UMAP1", "UMAP2", "NMF_Group")
  if (!is.data.frame(dat) || !all(required %in% names(dat))) {
    stop("Embedding '",nm,"' lacks required columns: ",paste(setdiff(required,names(dat)),collapse=", "))
  }
  if (anyNA(dat[[CONFIG$PATIENT_ID_COL]]) || anyDuplicated(dat[[CONFIG$PATIENT_ID_COL]])) {
    stop("Missing or duplicate source patient identifiers in ", nm)
  }
}
explorer_log("Loaded existing Script 16 results; no NMF/UMAP will be recalculated")
explorer_log("Raw rows: ", nrow(source_obj$embeddings$raw),
             "; normalized rows: ", nrow(source_obj$embeddings$normalized))

# ==============================================================================
# SECTION 14.04 — FIELD TYPES AND PALETTES
# ==============================================================================
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
# SECTION 14.05 — INTERACTIVE HTML ASSEMBLY AND SAVE
# ==============================================================================
make_unified_explorer <- function(embeddings) {
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
  approved_colors <- unique(as.character(approved_colors))
  assert_safe_fields(approved_colors, "color")
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
  assert_safe_fields(approved_hover, "tooltip")
  missing_hover <- setdiff(approved_hover, names(raw))
  if (length(missing_hover)) explorer_log("Unavailable hover columns: ", paste(missing_hover, collapse = ", "), level = "WARNING")
  approved_hover <- intersect(approved_hover, names(raw))
  # Only needed columns are serialized into the public-facing HTML document.
  data_cols <- unique(c("NMF_Group", available_colors, approved_hover))
  assert_safe_fields(data_cols, "embedded")
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

  # Check the actual serialized field names. Do NOT add USI, survival, or identifying
  # variables to HTML even if those columns exist in the source RDS.
  assert_safe_fields(names(export), "FINAL HTML")
  if (anyNA(export$NMF_Group) || any(!nzchar(export$NMF_Group))) {
    stop("NMF_Group is missing for one or more patients")
  }
  # Random display labels are independent of the source USI values. The
  # translation is held in memory only and is NEVER exported to HTML or disk.
  # These names are pseudonyms, NOT a guarantee of anonymization.
  set.seed(EXPLORER_CONFIG$PUBLIC_ID_SEED)
  export$Public_ID <- sample(sprintf("AML_%04d", seq_len(nrow(export))), nrow(export))

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
  destination <- file.path(EXPLORER_DIR, "Unified_UMAP_Explorer_PUBLIC.html")
  # A truly stand-alone Github Pages page should not require a separate
  # dependencies directory. Raise an error rather than silently publishing a
  # broken page if Pandoc/selfcontained conversion is unavailable.
  htmlwidgets::saveWidget(widget, destination,
                          selfcontained = isTRUE(EXPLORER_CONFIG$SELF_CONTAINED),
                          title = "Pediatric AML NMF — Interactive UMAP Explorer")
  explorer_record_output(destination)
  explorer_log("Unified explorer created: ", destination)
  explorer_log("Embedded PSEUDONYMIZED rows: ", nrow(export))
  explorer_log("Public color variables: ", paste(available_colors, collapse = ", "))
  explorer_log("Public tooltip-only variables: ", paste(setdiff(approved_hover, available_colors), collapse = ", "))
  explorer_log("IMPORTANT: Review the HTML and embedded patient-level attributes ",
           "before public GitHub publication.", level = "WARNING")
  # Private audit contains field names and counts only, not individual values.
  audit <- data.frame(field = names(export),
                      storage_type = vapply(export, function(z) class(z)[1L], character(1)),
                      missing_n = vapply(export, function(z) sum(is.na(z)), integer(1)),
                      stringsAsFactors = FALSE)
  utils::write.csv(audit, file.path(EXPLORER_DIR, "Public_Explorer_Field_Audit_PRIVATE.csv"), row.names = FALSE)
  explorer_log("Field audit saved privately alongside the candidate HTML")
  invisible(destination)
}


# ==============================================================================
# SECTION 14.06 — EXECUTE AND VERIFY DELIVERABLE
# ==============================================================================
# The only output intended for GitHub is the reviewed *_PUBLIC.html file.
# The field audit and execution log remain in your private output directory.
export_error <- tryCatch({
  output_file <- make_unified_explorer(source_obj$embeddings)
  explorer_log("Completed public-candidate explorer export")
  explorer_log("Review all embedded attributes before uploading to GitHub Pages",
               level="WARNING")
  NULL
}, error = function(e) {
  explorer_log("ERROR: ",conditionMessage(e),level="ERROR")
  e
})
if (inherits(export_error,"error")) stop(export_error)
# ==============================================================================
# END — STANDALONE PUBLIC EXPLORER EXPORTER
# ==============================================================================
