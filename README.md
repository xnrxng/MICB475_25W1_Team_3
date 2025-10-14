# [Project title]

*Authors: Rui Xiang Yu, Houria Afshar Moghaddam, Brooke Macleod, Quinlan Torstensen*

This repository houses the scripts, reports, and results from the project [name of project]. This project's corresponding paper is [link to paper].

The objective of this project is to [fill]. More information on the contents of this repository is below.

## Timeline and results

[gantt chart]

The project's report/digital lab notebook can be found in the `reports/` folder. Available formats are Quarto MarkDown [link], PDF [link], and HTML [link].

## Team meetings

Here are the notes and the agenda for every team meeting:

| Month | Week 1 | Week 2 | Week 3 | Week 4 |
|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
| October | [October 1st](https://github.com/xnrxng/MICB475_25W1_Team_3/blob/main/team_meetings/oct_1st_2025_meeting.md%3E) | [October 8th](https://github.com/xnrxng/MICB475_25W1_Team_3/blob/main/team_meetings/oct_8th_2025_meeting.md) | [October 15th](https://github.com/xnrxng/MICB475_25W1_Team_3/blob/main/team_meetings/oct_15th_2025_meeting.md) | link |
| November | link | link | link | link |

## Reproducibility

To reproduce the results in this repository, you should have the following set up:

1.  Git clone the repository:

    ``` bash
    git clone https://github.com/xnrxng/MICB475_25W1_Team_3.git
    ```

2.  `cd` to the root of the repository.

3.  Create the conda environment:

    ``` bash
    conda env create -f environment.yml
    ```

4.  Have the FASTQ files downloaded (need to add more info to this).

5.  Have the manifest.tsv file ready (add more info).

After this is set-up, you can re-run ALL the code using `make all` . Please note that this might be a time-consuming and memory-intensive process. We recommend starting a double-pane `tmux` session where `htop` can be run, while `make all` is running.
