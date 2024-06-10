FROM jupyter/minimal-notebook:lab-3.6.3

LABEL maintainer "Eppelt, Jonas jonas.eppelt@student.kit.edu"

USER root

### run apt-get installs
RUN apt-get update -y && apt-get upgrade -y && \
     apt-get install -y nano vim pandoc dvipng git texlive-xetex texlive texlive-lang-german texlive-latex-extra openssh-client wget sudo&& \
     apt-get clean && rm -rf /var/lib/apt/lists/*

### upgrading some packages
RUN python3.10 -m pip install metakernel ipython ipywidgets nbconvert notebook jupyterlab ipympl nbdime jupyterlab-git tqdm --upgrade --no-cache-dir


### copy extra start.sh for better mounting
COPY start.sh /usr/local/bin/


RUN conda install -c conda-forge ipympl ipywidgets matplotlib pandas pytables jupyterlab_rise
### switching user back
USER $NB_USER
