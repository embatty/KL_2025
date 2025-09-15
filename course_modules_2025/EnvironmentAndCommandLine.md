# Computational Practical: Setting up the environment for bioinformatics
## Table of Contents
1. [Introduction to bioinformatics software](#intro)
2. [Installing Conda](#conda)
3. [Download the course data files](#downloads)
4. [Set up Conda environments](#environments)

---

## Introduction to bioinformatics software <a name="intro"></a>
To use bioinformatic software, there are usually two options: running software tools on the command line, or using a graphic interface. Not all tools have a graphical interface set up, and the command line can be a powerful way to run bioinformatics software, so for this course we will mainly run software on the command line.

Most bioinformatics software will not run directly on a Windows computer. Instead, we need to set up a Unix command line where the software can be installed. The easiest way to set this up is to install [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) using the default settings.

For users of Mac computers, there is a built-in command line tool which can be found in the Applications/Utilities folder (open the "Terminal" program). However, if you have a modern Mac which runs on the Apple Silicon architecture, some software is not available. In this case, it can be helpful to [set up a dedicated Terminal program using emulation of older Macs](https://www.courier.com/blog/tips-and-tricks-to-setup-your-apple-m1-for-development) and use this to run your software.

## Installing Conda <a name="conda"></a>

To install bioinformatics software, we will use a tool call `conda`. Conda is an open-source environment and package manager, which will simplify the task of installing the software we need, and keeping track of software versions. Software packages can be found in conda channels, including a specific channel called 'bioconda' which includes many bioinformatics packages.

We will install conda using the miniforge installer, which allows us to install conda from the command line. To install conda, make sure you have a Unix terminal open (for Macs, this is the Terminal program, for Windows make sure you have Windows Subsystem for Linux open and not the Powershell terminal). Run the following command:
```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

This will download the correct installer. Then you can run the installer with this command line:

`bash Miniforge3-$(uname)-$(uname -m).sh`

## Download the course data files <a name="downloads"></a>

Download the course data files from [link]() and put them in a folder that you can see in your command line terminal.
(instructions for windows go here).

For each module, there is a folder which contains the data you will need to run the practical exercises for each module. There is also a subdirectory named `outputs` which has the output files produced in each module - if you have any trouble running any of the steps of the practical, you can look at the output anyway.

## Setting up conda environments <a name="environments"></a>

A conda environment is a convenient way to install a set of software tools so that you can always run the same software versions, and you can also send a configuration file to any other user, ensuring that your research is reproducible.

We have provided conda environment files, which will automatically set up a conda environment for you with all the software you need for each module.

All the environment files are in the `EnvironmentFile` subdirectory. For each module there is one environment file, which is named after the module and ends with the `.yaml` file extension.



```bash
docker run -it --mount type=bind,source=C:\Users\4M\Desktop\data,target=/home/data amr:Dockerfile
```

docker run -it --mount type=bind,src=/Users/Elizabeth/dockertest,dst=/home
Navigate to directory `home/data/`, download the course data, and create a new directory for this practical named `cp6`:

```bash
cd home/data/
git clone https://github.com/WCSCourses/AMR_2025/
mkdir cp6
```
