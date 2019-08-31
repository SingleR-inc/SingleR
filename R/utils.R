#' @importFrom SummarizedExperiment assay
#' @importFrom DelayedMatrixStats rowAnyNAs
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom DelayedArray DelayedArray
.to_clean_matrix <- function(x, assay.type, check.missing, msg="x") {
    if (is.null(rownames(x))) {
        stop(sprintf("'%s' must have row names", msg))
    }
    if (is(x, "SummarizedExperiment")) {
        x <- assay(x, i=assay.type)
    }

    # Stripping out genes with NA's from 'x'.
    if (check.missing) {
        discard <- rowAnyNAs(DelayedArray(x))
        if (any(discard)) {
            warning(sprintf("'%s' contains rows with missing values", msg))
            x <- x[!discard,,drop=FALSE]
        }
    }

    x
}

.make_heatmap_annotation_colors <- function(args, show.pruned) {
    # Create pheatmap annotations_colors dataframe
        # list of character vectors, all named.
            # vector names = annotation titles
            # vector members' (colors') names = annotation identities 
    
    # Extract a default color-set
    annotation.colors.d <- .make_heatmap_colors_discrete(show.pruned)
    annotation.colors.c <- .make_heatmap_colors_continuous()
    
    # Initiate variables
    next.color.index.discrete <- 1
    next.color.index.numeric <- 1
    
    col_colors <- NULL
    row_colors <- NULL
    
    # Columns First (if there)
    if (!is.null(args$annotation_col)) {

        for (i in seq_len(ncol(args$annotation_col))){
            
            # Determine the distinct contents of the first annotation
            in.this.annot <- levels(as.factor(args$annotation_col[,i]))
            
            # Make new colors
            if(!is.numeric(args$annotation_col[,i])){
                # Take colors for each, and name them.
                new.colors <- annotation.colors.d[
                    seq_along(in.this.annot) + next.color.index.discrete - 1
                    ]
                names(new.colors) <- in.this.annot
                
                next.color.index.discrete <-
                    next.color.index.discrete + length(in.this.annot)
            } else {
                # Make a 100 color split as in pheatmap code.
                a <- cut(args$annotation_col[,i], breaks = 100)
                # Assign to colors.
                this.ramp <- annotation.colors.c[next.color.index.numeric]
                new.colors <-
                    grDevices::colorRampPalette(c("white",this.ramp))(100)[a]
                
                next.color.index.numeric <- next.color.index.numeric + 1
            }
            # Add the new colors as the list
            col_colors <- c(
                col_colors,
                list(new.colors))
        }
        names(col_colors) <- names(args$annotation_col)
    }
    
    # Rows Second (if there)
    if (!is.null(args$annotation_row)) {

        for (i in seq_len(ncol(args$annotation_row))){
            
            # Determine the distinct contents of the first annotation
            in.this.annot <- levels(as.factor(args$annotation_row[,i]))
            
            # Make new colors
            if(!is.numeric(args$annotation_row[,i])){
                # Take colors for each, and name them.
                new.colors <- annotation.colors.d[
                    seq_along(in.this.annot) + next.color.index.discrete - 1
                    ]
                names(new.colors) <- in.this.annot
                
                next.color.index.discrete <-
                    next.color.index.discrete + length(in.this.annot)
            } else {
                # Make a 100 color split as in pheatmap code.
                a <- cut(args$annotation_row[,i], breaks = 100)
                # Assign to colors.
                this.ramp <- annotation.colors.c[next.color.index.numeric]
                new.colors <-
                    grDevices::colorRampPalette(c("white",this.ramp))(100)[a]
                
                next.color.index.numeric <- next.color.index.numeric + 1
            }
            # Add the new colors as the list
            row_colors <- c(
                row_colors,
                list(new.colors))
        }
        names(col_colors) <- names(args$annotation_col)
    }
    
    args$annotation_colors <- c(col_colors, row_colors)
    args
}

.make_heatmap_colors_discrete <- function(show.pruned) {
    # Creates a default vector of colors with 24*10 (overkill) options.
    annotation.colors <- rep(
        c(  # DittoSeq-v0.2.10 Colors
            "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
            "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
            "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
            "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C"),
        10)
    if (show.pruned) {
        annotation.colors <- c("white", annotation.colors)
    }
    annotation.colors
}

.make_heatmap_colors_continuous <- function() {
    # Creates a default vector of colors with 8*3 (overkill) options.
        # These represent max.colors for discrete color scales.
    annotation.colors <- rep(
        c(  # DittoSeq-v0.2.10 Colors, distinct order 
            "#B14380", "#A04700", "#005685", "#D5C711", "#007756",
            "#1C91D4", "#AD7700", "#4D4D4D", "#CC79A7", "#D55E00",
            "#0072B2", "#F0E442", "#009E73", "#56B4E9", "#E69F00",
            "#666666"),
        3)
    
    annotation.colors
}

