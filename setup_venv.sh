#!/bin/bash
# Script to setup the virtual environment and install dependencies

echo "Creating virtual environment..."
python3 -m venv shape-optimizer

echo "Activating environment and installing dependencies..."
source shape-optimizer/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

echo "Setup complete. Run 'source shape-optimizer/bin/activate' to start."