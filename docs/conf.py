project = "OpenMC Optimizer"
author = "OpenMC Optimizer contributors"
release = "0.1.0"

extensions = []
root_doc = "index"
templates_path = []
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_title = "OpenMC Optimizer Documentation"
html_theme_options = {
    "navigation_with_keys": True,
    "light_css_variables": {
        "color-brand-primary": "#075b70",
        "color-brand-content": "#087b96",
    },
    "dark_css_variables": {
        "color-brand-primary": "#61c5d8",
        "color-brand-content": "#8bd8e7",
    },
}
html_static_path = []
