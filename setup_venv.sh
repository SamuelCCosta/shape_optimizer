#!/bin/bash
# Script to setup the virtual environment and install dependencies

echo "Creating virtual environment..."
python3 -m venv venv

echo "Activating environment and installing dependencies..."
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

echo "Setup complete. Run 'source venv/bin/activate' to start."