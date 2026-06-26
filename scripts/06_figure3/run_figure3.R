# Figure 3 workflow driver.
# Run this script from the project root to rebuild the Figure 3 data tables,
# statistical summaries, and PNG panels in dependency order.
source(file.path("scripts", "06_figure3", "01_build_figure3_data.R"))
source(file.path("scripts", "06_figure3", "02_build_figure3_tables.R"))
source(file.path("scripts", "06_figure3", "03_plot_figure3.R"))
