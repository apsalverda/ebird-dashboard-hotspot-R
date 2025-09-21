# helper functions for ebird-dashbird-hotspot.qmd

read_ebird_data = function(my_filename = "MyEBirdData.csv"){
  read_csv(
    my_filename,
    col_types = cols(
      `Submission ID` = col_character(),
      `Common Name` = col_character(),
      `Scientific Name` = col_character(),
      `Taxonomic Order` = col_double(),
      Count = col_character(), #!
      `State/Province` = col_character(),
      County = col_character(),
      `Location ID` = col_character(),
      Location = col_character(),
      Latitude = col_double(),
      Longitude = col_double(),
      Date = col_date(format = ""),
      Time = col_time(format = ""),
      Protocol = col_character(),
      `Duration (Min)` = col_double(),
      `All Obs Reported` = col_double(),
      `Distance Traveled (km)` = col_double(),
      `Area Covered (ha)` = col_logical(),
      `Number of Observers` = col_double(),
      `Breeding Code` = col_character(),
      `Observation Details` = col_character(),
      `Checklist Comments` = col_character(),
      `ML Catalog Numbers` = col_character()
    ),
    na = c("")
  ) %>%
    janitor::clean_names() %>%
    select(submission_id, date, everything()) %>%
    arrange(date, time) %>%
    mutate(
      year = factor(year(date)),
      month = month(date),
      day_of_year = yday(date),
      us_season = case_when(
        month < 3 | month == 12~ "winter",
        month < 6 ~ "spring",
        month < 9 ~ "summer",
        month < 12 ~ "fall",
        TRUE ~ NA
      ),
      us_season = factor(us_season, levels = c("spring", "summer", "fall", "winter")),
      month = factor(month(date, label = TRUE)),
      # convert "X" to -1 so count can be stored as an integer
      count = as.numeric(ifelse(count == "X", "-1", count))
    ) |>
    filter(!str_detect(common_name, "sp\\."))
}

ordinal_date_suffix = Vectorize(
  function(date_str, year = TRUE, markdown_subscript = FALSE) {
    date_obj = as.Date(date_str, "%Y-%m-%d")
    day_suffix = c("st", "nd", "rd", rep("th", 17), "st", "nd", "rd", rep("th", 7), "st")
    if (markdown_subscript == TRUE){ day_suffix = paste0("^", day_suffix, "^") }
    formatted_date = paste0(
      month(date_obj, label = TRUE, abbr = FALSE),
      " ",
      day(date_obj),
      day_suffix[day(date_obj)],
      ifelse(year == TRUE, paste0(
        ", ",
        year(date_obj)
      ),
      ""
      )
    )
    return(formatted_date)
  }
)

my_percent = function(proportion){
  percentage = proportion * 100
  formatted_percentage = sprintf("%05.1f%%", percentage)
  return(formatted_percentage)
}

bg = function(start, end, color, ...) {
  paste(
    "linear-gradient(90deg,transparent ",
    my_percent(start),
    ", ",
    color,
    my_percent(start),
    ", ",
    color,
    my_percent(end),
    ", ",
    light_gray,
    my_percent(end),
    ")"
  )
}

xnormalize = function(x, y){
  y + (x / 100) * (1 - y)
}

my_color_bar = function(...){
  formatter(
    "span",
    style = function(x) style(
      display = "inline-block",
      "text-align" = "left",
      "width" = "100%",
      "background" = bg(0, xnormalize(x, 0), ebird_green)
    )
  )
}

year_list_ecdf = function(my_data, year_from = 2025, year_to = 2025){
  first_letter_months = function(x) {
    month_num = as.numeric(format(x, "%m"))
    first_letters = substr(month.abb, 1, 1)
    return(first_letters[month_num])
  }
  cum_counts = function(year_dat){
    year_dat %>%
      filter(year == year_from | year == year_to) %>%
      group_by(year, common_name) %>%
      slice(1) %>%
      arrange(year, date) %>%
      group_by(year) %>%
      mutate(
        total = row_number(),
        date = as.Date(str_replace(as.character(date), as.character(year_from), as.character(year_to)))
      ) %>%
      ungroup() |>
      select(date, total, year)
  }
  my_dat_from = cum_counts(my_data |> filter(year == year_from))
  my_dat_to = cum_counts(my_data |> filter(year == year_to))

  # add data point for last day of year
  if (!(month(max(my_dat_from$date)) == 12 & month(max(my_dat_from$date)) == 31)){
    my_dat_from = bind_rows(
      my_dat_from,
      my_dat_from |>
        last() |>
        mutate(date = ceiling_date(date, "year") - 1)
    )
  }
  my_plot =
    my_dat_from |>
    ggplot(aes(x = date, y = total))
  if (year_from != year_to){
    my_plot =
      my_plot +
      geom_step(
        color = "lightgray",
        linewidth = .3
      ) +
      geom_ribbon(
        aes(ymin = 0, ymax = total),
        fill = "lightgray",
        alpha = .2
      )
  }
  my_plot =
    my_plot +
    geom_step(
      data = my_dat_to,
      aes(x = date, y = total),
      color = "#2780e3",
      linewidth = .5
    ) +
    scale_linewidth_manual(values = c(.5, 1)) +
    scale_x_date(
      limits = c(as.Date(paste0(year_to, "-01-01")), as.Date(paste0(year_to, "-12-31"))),
      labels = first_letter_months,
      #      date_labels = c("J", substr(month.abb, 1, 1)), # HACK
      date_breaks = "1 month",
      expand = expansion(mult = c(.005, .005))
    ) +
    scale_y_continuous(expand = expansion(mult = c(.01, .05))) +
    labs(
      title = "Number of species reported",
      x = "",
      y = ""
    ) +
#    theme_gray(base_size = 14) +
    theme(
      text = element_text(size = 14),
      panel.background = element_blank(),
      # legend.position = "right",
      # legend.position.inside = c(.9, .12),
      # legend.background = element_rect(color = NA, fill = NA),
      # legend.box.background = element_rect(color = NA, fill = NA),
      # legend.title = element_blank(),
      plot.margin = margin(0, 10, 0, 0, "pt")
    )
  return(my_plot)
}

my_rects = function(my_plot){
  # generated vertical bands for seasons, for the Checklist plot
  rect_colors = c("white", "lightgray")
  alpha_level = .1
  rect_tibble = tibble()
  ls = layer_scales(my_plot)
  min_x = as.Date(ls$x$get_limits()[1])
  max_x = as.Date(ls$x$get_limits()[2])
  min_y = ls$y$get_limits()[1]
  max_y = ls$y$get_limits()[2]
  current_date = min_x
  rect_start_date = min_x
  rect_number = 1
  while(current_date < max_x){
    current_date = current_date + 1
    if (month(current_date) %in% c(3, 6, 9, 12) & day(current_date) == 1){
      rect_tibble =
        bind_rows(
          rect_tibble,
          tibble(
            xmin = rect_start_date,
            xmax = current_date,
            ymin = min_y,
            ymax = max_y,
            fill = rect_colors[1 + rect_number %% 2],
            alpha = alpha_level
          )
        )
      rect_start_date = current_date
      rect_number = rect_number + 1
    }
  }
  rect_tibble =
    bind_rows(
      rect_tibble,
      tibble(
        xmin = rect_start_date,
        xmax = current_date,
        ymin = min_y,
        ymax = max_y,
        fill = rect_colors[1 + rect_number %% 2],
        alpha = alpha_level
      )
    )
  return(rect_tibble)
}
