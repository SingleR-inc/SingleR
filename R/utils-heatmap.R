.make_heatmap_annotation_colors <- function(args, show.pruned) {
    # Create pheatmap annotations_colors dataframe
        # list of character vectors, all named.
            # vector names = annotation titles
            # vector members' (colors') names = annotation identities 
    
    # Extract a default color-set
    annotation.colors.d <- .make_heatmap_colors_discrete(show.pruned)
    annotation.colors.n <- .make_heatmap_colors_numeric()
    
    # Initiate variables
    next.color.index.discrete <- 1
    next.color.index.numeric <- 1
    col_colors <- NULL
    row_colors <- NULL
    
    # Columns First (if there)
    if (!is.null(args$annotation_col)) {
        dfcolors_out <- .pick_colors_for_df(
            args$annotation_col,
            next.color.index.discrete, next.color.index.numeric,
            annotation.colors.d, annotation.colors.n)
        col_colors <- dfcolors_out$df_colors
        next.color.index.discrete <- dfcolors_out$next.color.index.discrete
        next.color.index.numeric <- dfcolors_out$next.color.index.numeric
    }
    
    # Rows Second (if there)
    if (!is.null(args$annotation_col)) {
        dfcolors_out <- .pick_colors_for_df(
            args$annotation_col,
            next.color.index.discrete, next.color.index.numeric,
            annotation.colors.d, annotation.colors.n)
        row_colors <- dfcolors_out$df_colors
        next.color.index.discrete <- dfcolors_out$next.color.index.discrete
        next.color.index.numeric <- dfcolors_out$next.color.index.numeric
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

.make_heatmap_colors_numeric <- function() {
    # Creates a default vector of colors with 8*3 (overkill) options.
        # These represent max.colors for discrete color scales.
    rep(
        c(  # DittoSeq-v0.2.10 Colors, distinct order 
            "#B14380", "#A04700", "#005685", "#D5C711", "#007756",
            "#1C91D4", "#AD7700", "#4D4D4D", "#CC79A7", "#D55E00",
            "#0072B2", "#F0E442", "#009E73", "#56B4E9", "#E69F00",
            "#666666"),
        3)
}

# Interpret annotations dataframe,
# Pick, name, and add colors.
.pick_colors_for_df <- function(
    annotation_df,
    next.color.index.discrete, next.color.index.numeric,
    annotation.colors.d, annotation.colors.n
    ) {
    df_colors <- NULL
    for (i in seq_len(ncol(annotation_df))){
        
        # Determine the distinct contents of the first annotation
        in.this.annot <- levels(as.factor(annotation_df[,i]))
        
        # Make new colors
        if(!is.numeric(annotation_df[,i])){
            # Take colors for each, and name them.
            new.colors <- annotation.colors.d[
                seq_along(in.this.annot) + next.color.index.discrete - 1
                ]
            names(new.colors) <- in.this.annot
            
            next.color.index.discrete <-
                next.color.index.discrete + length(in.this.annot)
        } else {
            # Make a 100 color split as in pheatmap code.
            a <- cut(
                annotation_df[order(annotation_df[,i]),i],
                breaks = 100)
            # Assign to colors.
            this.ramp <- annotation.colors.n[next.color.index.numeric]
            new.colors <-
                grDevices::colorRampPalette(c("white",this.ramp))(100)[a]
            
            next.color.index.numeric <- next.color.index.numeric + 1
        }
        # Add the new colors as the list
        df_colors <- c(
            df_colors,
            list(new.colors))
    }
    names(df_colors) <- names(annotation_df)
    list(df_colors = df_colors,
         next.color.index.discrete = next.color.index.discrete,
         next.color.index.numeric = next.color.index.numeric)
}


