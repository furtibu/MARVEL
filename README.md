# What is MARVEL?

The acronym MARVEL stands for Measured Active Rotational-Vibrational Energy Levels. MARVEL is based on the concept of spectroscopic networks. The code MARVEL is written in C++ and it is designed for 
critically evaluating and validating experimental transition wavenumbers and uncertainties collected from the literature as well as for inverting the wavenumber information in order to obtain the best 
possible energy levels with attached uncertainties in a highly efficient way. MARVEL simultaneously processes all the available assigned experimental lines and determines the associated energy levels for 
the chosen molecule. In a way MARVEL can invert a 100,000 by 100,000 matrix in a fraction of a second. While the energy levels obtained could be called empirical, in order to emphasize the process they 
were derived from the evaluated and validated energy levels are usually called MARVEL energy levels.

## What is MARVEL?

| Molecule         | Isotopologue             | \#lines  | \#levels | Year | Reference                                |
|------------------|--------------------------|----------|----------|------|------------------------------------------|
| AlH         | <sup>27</sup>AlH               | 682      | 259      | 2018 | Yurchenko \ea\     \cite{AlH            |
| BeH         | <sup>9</sup>BeH                | 2201     | 1264     | 2018 | Darby-Lewis \ea\   \cite{BeH            |
|                  | <sup>9</sup>BeD                | 2605     | 1495     | 2018 | Darby-Lewis \ea\   \cite{BeH            |
|                  | <sup>9</sup>BeT                | 538      | 215      | 2018 | Darby-Lewis \ea\   \cite{BeH            |
| C<sub>2</sub>       | <sup>12</sup>C<sub>2</sub>              | 31\,323  | 7047     | 2020 | McKemmish \ea\     \cite{McKemmish_2020 |
| CH          | <sup>12</sup>CH                | 6348     | 1521     | 2022 | Furtenbacher \ea\  \cite{22FuHeTeCsa    |
| CN          | <sup>12</sup>C<sup>14</sup>N         | 40\,516  | 8083     | 2020 | Syme \ea           \cite{Syme_2020      |
| CP          | <sup>12</sup>C<sup>31</sup>P         | 3264     | 948      | 2021 | Qin \ea\           \cite{Qin_2021       |
| NH          | <sup>14</sup>NH                | 3002     | 1058     | 2019 | Darby-Lewis \ea\   \cite{19DaShJoKh     |
| NO          | <sup>14</sup>N<sup>16</sup>O         | 11\,136  | 4106     | 2017 | Wong \ea\          \cite{Wong_2017      |
| O<sub>2</sub>        | <sup>16</sup>O<sub>2</sub>              | 30\,671  | 15\,946  | 2019 | Furtenbacher \ea\  \cite{19FuHoKoSo     |
| OH          | <sup>16</sup>OH                | 15\,938  | 1624     | 2022 | Furtenbacher \ea\  \cite{22FuHeTeCsa    |
| PN          | <sup>31</sup>P<sup>14</sup>N         | 1735     | 61       | 2021 | Qin \ea\           \cite{Qin_2021_PN    |
| TiO         | <sup>48</sup>Ti<sup>16</sup>O        | 56\,240  | 10\,761  | 2019 | McKemmish \ea\     \cite{McKemmish_2019 |
| ZrO         | <sup>90</sup>Zr<sup>16</sup>O        | 22\,549  | 8088     | 2018 | McKemmish \ea\     \cite{18McBoGoSh     |
| CaOH        | <sup>40</sup>Ca<sup>16</sup>OH       | 3~240    | 1~955    | 2020 | Wang \ea\          \cite{Wang_2020      |
| H<sub>2</sub> O      | H<sub>2</sub> <sup>16</sup>O             | 289\,070 | 19\,204  | 2020 | Furtenbacher \ea\  \cite{20FuToTePo     |
|                  | H<sub>2</sub> <sup>17</sup>O             | 27\,045  | 5278     | 2020 | Furtenbacher \ea\  \cite{20FuCoTeYu_2   |
|                  | H<sub>2</sub> <sup>18</sup>O             | 66\,166  | 6865     | 2020 | Furtenbacher \ea\  \cite{20FuCoTeYu_2   |
|                  | HD<sup>16</sup>O               | 54\,740  | 8819     | 2010 | Tennyson \ea\      \cite{10TeBeBrCa     |
|                  | HD<sup>17</sup>O               | 485      | 162      | 2010 | Tennyson \ea\      \cite{10TeBeBrCa     |
|                  | HD<sup>18</sup>O               | 8728     | 1864     | 2010 | Tennyson \ea\      \cite{10TeBeBrCa     |
|                  | D<sub>2</sub> <sup>16</sup>O             | 53\,534  | 12\,269  | 2014 | Tennyson \ea\      \cite{14TeBeBrCa     |
|                  | D<sub>2</sub> <sup>17<sup>O             | 600      | 338      | 2014 | Tennyson \ea\      \cite{14TeBeBrCa     |
|                  | D<sub>2</sub> <sup>18</sup>O             | 12\,146  | 3350     | 2014 | Tennyson \ea\      \cite{14TeBeBrCa     |
| H<sub>2</sub> S      | H<sub>2</sub> <sup>32</sup>S             | 44\,325  | 7436     | 2018 | Chubb \ea\         \cite{18ChNaKeBa     |
| H<sub>3</sub> <sup>+</sup>   | H<sub>3</sub> <sup>+<sup>                | 1610     | 652      | 2013 | Furtenbacher \ea\  \cite{13FuSzMaFa     |
|                  | H<sub>2</sub> D<sup>+</sup>             | 195      | 86       | 2013 | Furtenbacher \ea\  \cite{13FuSzMaFa     |
|                  | D<sub>2</sub> H<sup>+</sup>             | 154      | 72       | 2013 | Furtenbacher \ea\  \cite{13FuSzMaFa     |
| HOCl        | H<sup>16</sup>O<sup>35</sup>Cl       | 20\,349  | 5760     | 2023 | Ecseri \ea\        \cite{23EcSiFuRa     |
|                  | H<sup>16</sup>O<sup>37</sup>Cl       | 10\,266  | 3933     | 2023 | Ecseri \ea\        \cite{23EcSiFuRa     |
| SO<sub>2</sub>       | <sup>32<sup>S<sup>16</sup>O<sub>2</sub>      | 40\,325  | 15\,130  | 2018 | T\'obi\'as \ea\    \cite{18ToFuCsNa     |
|                  | <sup>33</sup>S<sup>16</sup>O<sub>2</sub>      | 15\,647  | 5852     | 2018 | T\'obi\'as \ea\    \cite{18ToFuCsNa     |
|                  | <sup>34</sup>S<sup>16</sup>O<sub>2</sub>      | 31\,088  | 10\,893  | 2018 | T\'obi\'as \ea\    \cite{18ToFuCsNa     |
|                  | <sup>36</sup>S<sup>16</sup>O<sub>2</sub>      | 31       | --       | 2018 | T\'obi\'as \ea\    \cite{18ToFuCsNa     |
| NH<sub>3</sub>       | <sup>14</sup>NH<sub>3</sub>             | 29\,450  | 4961     | 2015 | Al Derzi \ea\      \cite{15AlFuYuTe     |
| C<sub>2</sub>H<sub>2</sub>  | <sup>12</sup>C<sub>2</sub>H<sub>2</sub>        | 37\,813  | 11\,213  | 2018 | Chubb \ea\         \cite{18ChJoFrCh     |
| H<sub>2</sub>CO     | H<sub>2</sub><sup>12</sup>C<sup>16</sup>O     | 39\,662  | 5029     | 2021 | Al-Derzi \ea\      \cite{21DeTeYuMe     |
| H<sub>2</sub>C<sub>2</sub>O | H<sub>2</sub><sup>12</sup>C<sub>2</sub><sup>16</sup>O | 3879     | 353      | 2011 | F\'abri \ea\       \cite{11FaMaFuNe     |
