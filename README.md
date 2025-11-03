# Exploring the Interplay Between Fibre Intake, Exercise, and Gut Microbiota in Modulating Cardiovascular Health in a Westernization Context

*Authors: Rui Xiang Yu, Houria Afshar Moghaddam, Brooke Macleod, Quinlan Torstensen*

*Supervisors: Chad Poloni, Dr. Evelyn Sun, Dr. Avril Metcalfe-Roach*

This repository houses the scripts, reports, and results from the project "Exploring the Interplay Between Fibre Intake, Exercise, and Gut Microbiota in Modulating Cardiovascular Health in a Westernization Context". This project's corresponding paper will be published in UJEMI.

More information on the contents of this repository is below.

## Timeline and results

![](results/gantt_chart.png)

The project's report/digital lab notebook can be found in the `reports/` folder. Available formats are Quarto MarkDown [link], PDF [link], and HTML [link].

## Team meetings

Here are the notes and the agenda for every team meeting:

| Month | Week 1 | Week 2 | Week 3 | Week 4 | Week 5 |
|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|
| October | [October 1st](https://github.com/xnrxng/MICB475_25W1_Team_3/blob/main/team_meetings/oct_1st_2025_meeting.md) | [October 8th](https://github.com/xnrxng/MICB475_25W1_Team_3/blob/main/team_meetings/oct_8th_2025_meeting.md) | [October 15th](https://github.com/xnrxng/MICB475_25W1_Team_3/blob/main/team_meetings/oct_15th_2025_meeting.md) | [October 22nd](https://github.com/xnrxng/MICB475_25W1_Team_3/blob/main/team_meetings/oct_22nd_2025_meeting.md) | [October 29th](https://github.com/xnrxng/MICB475_25W1_Team_3/blob/main/team_meetings/oct_29th_2025_meeting.md) |
| November | [November 5th](https://github.com/xnrxng/MICB475_25W1_Team_3/blob/main/team_meetings/nov_5th_2025_meeting.md) | link | link | link | link |

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

4.  Activate the conda environment:

    ``` bash
    conda activate qiime2-amplicon-2025.4
    ```

5.  Have the FASTQ files downloaded (need to add more info to this).

6.  Have the manifest.tsv file ready (add more info).

After this is set-up, you can re-run ALL the code using `make all` . Please note that this might be a time-consuming and memory-intensive process. We recommend starting a double-pane `tmux` session where `htop` can be run, while `make all` is running.
