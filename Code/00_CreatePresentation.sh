#!/bin/bash
jupyter nbconvert 00_Presentation.ipynb --to slides --no-prompt --TagRemovePreprocessor.remove_input_tags={\"to_remove\"}