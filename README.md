qkcnotebooks

--- virtuelle Umgebung erstellen (Python 3.6)
virtualenv qkclab
--- virtuelle Umgebung aktivieren
cd qkclab/scripts
activate
--- virtuelle Umgebung für Jupyter Notebook einrichten
pip install ipykernel
python -m ipykernel install --user --name=qkclab
--- Jupyter öffnen
jupyter notebook

--- installing modules
--- activate

pip install matplotlib
pip install ipywidgets