
#calculate kest
calculate_gcross_exact = function(config, path){
  #set correct path
  if(!(getwd() == config$paths$root)){
    setwd(config$paths$root)
  }
  #read in files
  df = fread(path)
  #split based on the classifier label colun
  if(!is.null(config$variables$tissue_class_label)){
    df_list = split(df, df[[config$variables$tissue_class_label]])
  } else {
    df_list = list("Overall" = df)
  }
  #identify windows of compartments
  windows = lapply(df_list, function(x){
    convexhull.xy(x[[config$variables$x_value]], 
                  x[[config$variables$y_value]])
  })
  #point pattern objects
  pp_objs = lapply(seq(df_list), function(l){
    x = df_list[[l]]
    marks = x[,unique(unlist(config$variables$bivar_pairs)), with = FALSE]
    marks[] = lapply(marks, as.factor) #really is unneeded but YOLO
    ppp(x[[config$variables$x_value]], 
        x[[config$variables$y_value]], window = windows[[l]],
        marks = marks)
  })
  gcross_res = lapply(seq(df_list), function(l){
    p = pp_objs[[l]]
    #pairs
    res = lapply(config$variables$bivar_pairs, function(markers){
      #find cells to subset to
      mark_cells = marks(p) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), ~as.numeric(as.character(.x))), #forcing 1/0 for positive negative
                      id = 1:dplyr::n()) %>%
        dplyr::filter(get(markers[1]) != get(markers[2])) %>%
        dplyr::select(!!markers, id) #selects in correct order matching markers - anchor, counted
      counts = colSums(mark_cells[,markers, with = FALSE])
      if(any(counts < 2)){
        out = data.frame(r = eval(parse(text = config$variables$radii_range)))  %>%
          mutate(Anchor = markers[1],
                 Counted = markers[2], 
                 .before = 1)
      } else {
        cells = mark_cells %>%
          tidyr::pivot_longer(dplyr::all_of(markers), names_to = "marker", values_to = "positive") %>%
          dplyr::filter(positive == 1)
        p2 = subset(p, cells$id)
        marks(p2) = factor(cells$marker, levels = markers)
        obs = Gcross(p2, r=eval(parse(text = config$variables$radii_range)), correction = "none") %>%
          as.data.frame(check.names = FALSE) %>%
          mutate(Anchor = markers[1],
                 Counted = markers[2], 
                 .before = 1)
        N = npoints(p)
        A_idx = cells$id[cells$marker == markers[1]]
        B_idx = cells$id[cells$marker == markers[2]]
        nA = length(A_idx)
        K = length(B_idx)
        radii = obs$r
        counts_all = compute_neighbor_counts_opt(p$x, p$y, radii) 
        counts_A   = counts_all[A_idx, , drop = FALSE]
        zero_hit_probs = binom_prob(n = N - 1, k = counts_all, x = K)
        Gcross_exact = 1 - (1 / N) * colSums(zero_hit_probs)
        obs$raw_exact = Gcross_exact
        return(obs)
        #permutations
        # set.seed(333)
        # perms = parallel::mclapply(1:10000, function(rp){
        #   r_cells = sample(1:npoints(p), length(cells$id)) #random cells
        #   p2 = subset(p, r_cells)
        #   marks(p2) = factor(sample(cells$marker), levels = markers)
        #   perm = Gcross(p2, r=eval(parse(text = config$variables$radii_range)), correction = "raw") %>%
        #     as.data.frame(check.names = FALSE) %>%
        #     mutate(Anchor = markers[1],
        #            Counted = markers[2], 
        #            .before = 1)
        #   return(perm)
        # }, mc.cores = 16, mc.preschedule = FALSE) %>%
        #   do.call(dplyr::bind_rows, .) %>%
        #   group_by(Anchor, Counted, r) %>%
        #   summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
      }
      
    }) %>%
      do.call(bind_rows, .)
    
    if(!is.null(config$variables$tissue_class_label)){
      res %>%
        mutate(!!config$variables$tissue_class_label := names(df_list)[l]) %>%
        return()
    } else {
      return(res)
    }
    
  }) %>%
    do.call(bind_rows, .) %>%
    dplyr::mutate(!!config$variables$sample_id := unique(df[[config$variables$sample_id]]),
                  .before = 1)
  fwrite(gcross_res,
         file.path(config$paths$output, 'metrics/gcross_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".csv.gz")),
         compress = "gzip")
  
  return("Success")
  
}

#plot it
plot_gcross_exact = function(config, path){
  #set correct path
  setwd(config$paths$root)
  df = fread(path, data.table = FALSE) %>%
    mutate(r = as.numeric(r))
  if(!is.null(config$variables$tissue_class_label)){
    c_num = length(unique(df$Anchor))
    p = df %>%
      ggplot() + 
      geom_line(aes(x = r, y = raw - raw_exact, color = Counted)) +
      facet_grid(get(config$variables$tissue_class_label) ~ Anchor) +
      theme_classic()
    
    pdf(file.path(config$paths$output, 'figures/metrics/gcross_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 7, width = (5*c_num)+0.5)
    print(p)
    dev.off()
  } else {
    c_num = length(unique(df$Anchor))
    p = df %>%
      ggplot() + 
      geom_line(aes(x = r, y = rs - theo, color = Counted)) +
      facet_grid(. ~ Anchor) +
      theme_classic()
    
    pdf(file.path(config$paths$output, 'figures/metrics/gcross_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 7, width = (5*c_num)+0.5)
    print(p)
    dev.off()
  }
  
  return("Success")
}
