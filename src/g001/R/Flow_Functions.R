
#' Get overall experiment time
#'
#' @param exp_to_do experiment FCS file
#'
#' @return tibble with time and tube
#' @export
#'
get_exp_time <- function(exp_to_do) {
  temp_times <- exprs(gh_pop_get_data(experimental_samples_global[[exp_to_do]]))[,'Time']

  tibble(name = exp_to_do,
         Min_Time = min(temp_times),
         Max_Time = max(temp_times))
}



#' Getting Counts at specific times
#'
#' @param exp_to_do experiment FCS file
#' @param inx_to_do index FCS file
#'
#' @return counts for a specific time range
#' @export
#'
run_inx_stats_fun <- function(exp_to_do, inx_to_do) {
  inx_time <- exprs(gh_pop_get_data(index_samples_global[[inx_to_do]]))[,'Time']

  if (diff(range(inx_time)) == 0)
    cat('Warning: No Difference Between INX file', inx_to_do,
        'Min and Max Times', sep = '\n')


  gh_pop_get_stats_tfilter(experimental_samples_global[[exp_to_do]],
                           tfilter = list(min(inx_time), max(inx_time))) %>%
    mutate(Min_Time = min(inx_time),
           Max_Time = max(inx_time),
           INX_File = inx_to_do)
}



## Reading in the diva table (csv files) to compare to our counts/percents
read_diva_table <- function(f){
  input_data <-
    as_tibble(data.table::fread(f, fill = TRUE, sep = ","))
  split_at <- which(input_data$V1 == "Population")
  header <- input_data[1:(split_at - 1), ]
  body <- input_data[split_at:nrow(input_data), ]

  header <-
    header %>% dplyr::select(V1, V2) %>% tidyr::spread(V1, V2) %>%
    mutate(`Record Date` = as.Date(`Record Date`,
                                   tryFormats = c("%b %d, %Y", "%m/%d/%y"),
                                   optional = TRUE))

  colnames(body) <- body[1, ]
  body <- body[-c(1L), ]


  body <-
    body %>%
    mutate(
      count_diva = as.numeric(`#Events`),
      proportion_diva = as.numeric(`%Parent`),
    )


  # Some diva output csv lack the Parent Name column.
  # We handle that here.
  #
  if ("Parent Name" %in% colnames(body)) {
    body <-
      body %>%
      select(Population, Parent = `Parent Name`, count_diva, proportion_diva) %>%
      # tidy the input data so that Live/Dump is represented as Live|Dump
      mutate(Population = gsub("/","|",Population),
             Parent = gsub("/","|",Parent))

    body$Parent <- gsub("All Events", "root", body$Parent)
  } else {
    body <-
      body %>%
      select(Population, count_diva, proportion_diva) %>%
      # tidy the input data so that Live/Dump is represented as Live|Dump
      mutate(Population = gsub("/","|",Population))
  }

  output_data <- cbind(header, body) %>%
    # remove the row where Parent is empty and population is "All Events"
    dplyr::filter(Population != "All Events") %>%
    # Replace NA with 0, since it's due to 0 in the denominator when computing the proportion
    replace_na(replace = list(proportion_diva = 0)) %>%
    mutate("Tube Name" = gsub("V","",`Tube Name`))

  output_data
}


singel_gate_fix <- function(ws_here, gate_name, color_in = NA) {
  if (any(gs_get_pop_paths(ws_here, path = 'auto') == gate_name)) {
    tmp_gate <- try(gs_pop_get_gate(ws_here, gate_name))

    if (!inherits(tmp_gate, "try-error")) {
      # if color_in Na do both dim
      if (is.na(color_in))
        color_in = map_chr(tmp_gate[[1]]@parameters,~.x@parameters)


      for (i_name in names(tmp_gate)) {
        for (tmp_color in color_in) {
          # rectangle gate
          if (any(slotNames(tmp_gate[[i_name]]) == 'min')) {
            #Setting min to -10, which should be lower than any point
            if (!is.na(tmp_gate[[i_name]]@min[tmp_color]))
              tmp_gate[[i_name]]@min[tmp_color] = pmin(-10, tmp_gate[[i_name]]@min[tmp_color])
          }
          else {
            # polygonal gate
            tmp_boundaries <- tmp_gate[[i_name]]@boundaries[, tmp_color]
            tmp_index <- tmp_boundaries == min(tmp_boundaries)
            #Setting min to -10, which should be lower than any point
            tmp_gate[[i_name]]@boundaries[tmp_index, tmp_color] = pmin(rep(-10, sum(tmp_index)), min(tmp_boundaries))
          }
        }
        gh_pop_set_gate(ws_here[[i_name]], gate_name, tmp_gate[[i_name]])
      }
    }
    ws_here
  } else {
    cat("Warning: ", gate_name, "not present for manual adjusting\n")
    ws_here
  }
}


# Function to manually set boundaries to capture negative points because Diva messed up and incorrectly cut at 0.0
manually_adj_gates <- function(ws_in, IgD_color_in, visit_in) {
  gates_here <- gs_get_pop_paths(ws_in, path = 'auto')
  gates_here_full <- gs_get_pop_paths(ws_in, path = 'full')

  # IgD or IgDIgG gate
  if (any(gates_here == "IgD-")) {
    ws_in <- singel_gate_fix(ws_in, "IgD-")
  } else {
    ws_in <- singel_gate_fix(ws_in, "IgD-IgG+", color_in = IgD_color_in)
  }
  #Fixing eODKO11- gate
  if (any(gates_here == "eODKO11")) {
    ws_in <- singel_gate_fix(ws_in, "eODKO11")
  } else {
    ws_in <- singel_gate_fix(ws_in, "eODKO11-")
  }
  #Fixing GT8++KO- gate
  if (any(gates_here == "GT8++KO+")) {
    # VRC FNA samples mislabeled GT8++KO- as GT8++KO+
    ws_in <- singel_gate_fix(ws_in, "GT8++KO+")
  } else {
    ws_in <- singel_gate_fix(ws_in, "GT8++KO-")
  }
  #Live/Dump gate
  if (any(gates_here == "Live|Dump-"))
    ws_in <- singel_gate_fix(ws_in, "Live|Dump-")
  #Live/CD14-CD3-
  if (any(gates_here == "Live|CD14-CD3-"))
    ws_in <- singel_gate_fix(ws_in, "Live|CD14-CD3-")
  if (any(gates_here == "CD14-CD3-"))
    ws_in <- singel_gate_fix(ws_in, "CD14-CD3-")
  #CD14-CD3+ gate
  if (any(gates_here == "Live|CD14-CD3+"))
    ws_in <- singel_gate_fix(ws_in, "Live|CD14-CD3+", color_in = 'G660-A')
  if (any(gates_here == "CD14-CD3+"))
    ws_in <- singel_gate_fix(ws_in, "CD14-CD3+", color_in = 'G660-A')


  # Adding new gates for IgD- cells
  if (visit_in != 'V07A') {
    igd_igg_gate <- try(gs_pop_get_gate(ws_in, "IgD-IgG+"))
    if (class(igd_igg_gate[[1]]) == 'rectangleGate') {
      igd_cutoff <- map_dbl(igd_igg_gate,
                            ~.x@max[IgD_color_in]) %>%
        unique()
    } else {
    igd_cutoff <- map_dbl(igd_igg_gate,
                          ~.x@boundaries[,IgD_color_in] %>% max) %>%
      unique()
    }
    tmp_gate_dims <- list(c(-Inf, Inf),
                          c(-Inf, Inf))
    names(tmp_gate_dims) <- igd_igg_gate[[1]]@parameters %>%
      names()
    tmp_gate_dims[[IgD_color_in]][2] <- igd_cutoff

    rg <- rectangleGate(tmp_gate_dims, filterId = 'IgD-')
    igd_or_gc_nodeID <- gs_pop_add(ws_in,rg,parent = "CD19+CD20+")

    if (any(gates_here == "CD20hiCD38hi")) {
      gc_gate <- try(gs_pop_get_gate(ws_in, "CD20hiCD38hi"))
      igd_or_gc_nodeID <- gs_pop_add(ws_in,gc_gate,
                                    parent = 'IgD-')
    }

    if (any(gates_here == "GT8++noKO")) {
      gt8_nok0_gate <- try(gs_pop_get_gate(ws_in, "GT8++noKO"))
      gt8_nok0_nodeID <- gs_pop_add(ws_in,gt8_nok0_gate,
                                    parent = gs_get_pop_paths(ws_in)[igd_or_gc_nodeID])

      if (any(gates_here == "GT8++KO-")) {
        ep_specific_gate <- try(gs_pop_get_gate(ws_in, "GT8++KO-"))
        ep_specific_nodeIDs <- gs_pop_add(ws_in,ep_specific_gate,
                                          parent = gs_get_pop_paths(ws_in)[gt8_nok0_nodeID])
      }
    }
  }
  recompute(ws_in)
  ws_in
}


