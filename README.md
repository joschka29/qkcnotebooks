## qkcnotebooks

### virtuelle Umgebung erstellen (Python 3.6)
* pip install virtualenv
* virtualenv qkclab
### virtuelle Umgebung aktivieren
* cd qkclab/bin
* Linux: source activate
* Windows: activate
### virtuelle Umgebung für Jupyter Notebook einrichten (Kernel mit der venv erstellen)
* pip install ipykernel
* python -m ipykernel install --user --name=qkclab
* pip install matplotlib scipy numpy ipywidgets
### Jupyter öffnen und Kernel auswählen
* jupyter notebook
* Kernel > Change Kernel > qkclab

