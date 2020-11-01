@echo off

REM Run the main julia code
julia main.jl

REM Use the data from the julia run to create the figures
python figs.py

