# Globus

Globus is a file transfer utility used for research data management. It is supported by the Digital Research Aliance of Canada (Compute Canada) and UBC ARC clusters. There are two ways to transfer files: via the Globus web app, or via the Globus command line interface.

## Transfer via Globus Web App

The Globus web app can be accessed at http://globus.computecanada.ca (or simply https://globus.org). Follow these steps to transfer data:

0. Log in using your institutional login and grant the appropriate permissions

1. Once you have access to the Globus file manager, navigate to `Collections --> Search`. From here you can search for any Globus endpoint. For example, searching for `graham-dtn` (the data transfer node associated with the Graham cluster) will reveal the Globus endpoint `computecanada#graham-globus`. Select the desired endpoint and click `Open in File Manager`

**NOTE**: UBC ARC Sockeye has a collection called `ubcarc#sockeye`, but you should instead search for your specific collection (e.g. Matt Choptuik's Shared Chinook Allocation)

**NOTE**: You may want to bookmark your oft-used endpoints, which can be done from the file manager

2. Now in the file manager, click on the `Set Two Panels` button in the top right of the web app. A second panel should open up which will allow you to search for a second endpoint. For example, your could search for `cedar-dtn` if you wish to transfer between Graham and Cedar

3. Now transfer data by navigating to the appropriate directories and using the GUI to select the desired folders. You can also create and delete folders as desired. When you are ready to transfer, select 'Start' to begin transfer from the source to the destination

**NOTE**: You can also look at `Transfer & Timer Options` to adjust various options of the transfer

4. Once the transfer request has been submitted, you can monitor the transfer in the
`Activity` panel listed on the left of the page. At this point it is safe to close the web app. You should get an email when your transfer is complete

NOTE: If you don't want to deal with the Globus web app, you could also set up the Globus CLI on the ComputeCanada clusters (taking care to load the correct modules); see below

## Transfer via Globus CLI

Globus also has a CLI implemented in Python. Here I will use `venv` to install the CLI application. You can install `venv` with `sudo apt-get install python3-venv`, though it may already be installed by default

0. Create a virtual environment to install Globus CLI:
    ```
    > GLOBUSPATH=${HOME}/opt/globus-cli-virtualenv  # set as desired
    > python3 -m venv ${GLOBUSPATH}
    ```

    Activate the newly-created environment:
    ```
    > source ${GLOBUSPATH}/bin/activate
    ```

    Install Globus CLI into the virtual environment:
    ```
    > pip install globus-cli
    ```

**NOTE**: You will have to load the virtual environment every time you want to use the Globus CLI. If you don't like this, then just add the environment location to your path:
```
> export PATH=$PATH:${GLOBUSPATH}/bin
> echo 'export PATH=$PATH:${GLOBUSPATH}/bin'>>$HOME/.bashrc
```

**NOTE**: virtual environments are not easily movable due to how they set up their paths. This means that you should be sure of where you want to setup your virtual environment before you do it

1. Log on to Globus:
    ```
    > globus login --no-local-server
    ```
    Follow the instructions to authenticate the session

2. Search for globus endpoints:
    ```
    > globus endpoint search 'Choptuik'
    ID                                   | Owner                 | Display Name
    ------------------------------------ | --------------------- | -----------------------------------------
    <allocation ID will be shown here>   | <owner ID shown here> | Matt Choptuik's Chinook Allocation
    <allocation ID will be shown here>   | <owner ID shown here> | Matt Choptuik's Shared Chinook Allocation
    ```

    Alternatively, you can show your saved bookmarks:
    ```
    > globus bookmark list
    Name                                      | Bookmark ID              | Endpoint ID              | Endpoint Name                             | Path
    ----------------------------------------- | ------------------------ | ------------------------ | ----------------------------------------- | ------------
    computecanada#cedar-dtn                   | <bookmark ID shown here> | <allocation ID shown here> | computecanada#cedar-dtn                   | /home/mikin/
    Matt Choptuik's Shared Chinook Allocation | <bookmark ID shown here> | <allocation ID shown here> | Matt Choptuik's Shared Chinook Allocation | /mikin/
    ```

    Once you have identified the relevant endpoints, you probably want to save the ID into a shell variable:
    ```
    > ep1=<ID goes here>
    > ep2=<ID goes here>
    ```

3. You can perform various actions for each endpoint:

    List the files and directories::
    ```
    > globus ls ${ep1}:~/project
    mikin/
    test/
    ```

    Make a new directory:
    ```
    > globus mkdir ${ep1}:~/test_dir
    The directory was created successfully
    ```

    Rename a file or directory:
    ```
    > globus rename ${ep1}:~/test_dir ${ep1}:~/test_dir_renamed
    File or directory renamed successfully
    ```

    Transfer a single file from one endpoint to another:
    ```
    > globus transfer ${ep1}:~/file1.txt ${ep2}:~/file1_recv.txt --label "CLI single file"
    Message: The transfer has been accepted and a task has been created and queued for execution
    Task ID: 466a5962-dda0-11e6-9d11-22000a1e3b52
    ```

    Recursively transfer a directory from one endpoint to another:
    ```
    > globus transfer ${ep1}:~/test_dir ${ep2}:~/test_dir_recv --recursive --label "CLI single dir"
    Message: The transfer has been accepted and a task has been created and queued for execution
    Task ID: 47477b62-dda0-11e6-9d11-22000a1e3b52
    ```

    Use a .txt file to include multiple files in one transfer request:
    ```
    > cat in.txt
    file1.txt file1.txt # source path is followed by destination path
    file2.txt file2.txt # inline-comments are allowed in this file
    file3.txt file3.txt

    > globus transfer ${ep1}:~/test_dir ${ep2}:~/test_dir_recv --label "CLI Batch" --batch - < in.txt
    Message: The transfer has been accepted and a task has been created and queued for execution
    Task ID: 306900e0-dda1-11e6-9d11-22000a1e3b52
    ```
    
For further features of Globus CLI, see here:
https://docs.globus.org/cli/

