# SGGX-Convert-Tool
A tool which converts Zhao's distribution to SGGX

Zhao et al. represent fiber-like materials with angular Gaussian distributions D<sub>Zhao</sub> parameterized by a tangent direction ω = (x, y, z) and a roughness coefficient γ. This tool help user to generate SGGX data using the previous data in Zhao's format. 

Zhao's distribution require a precomputed table and some boost function and it has been realized in Mitsuba source code. Since it is somewhat laborious to extract these code due to complex dependencies, I made this tool.

This tool has some dependencies. To get these dependencies, install Mercurial first and input follow command in your command line.

<code>C:\>cd mitsuba</code>

<code>C:\mitsuba\>hg clone https://www.mitsuba-renderer.org/hg/dependencies_windows</code>

<code>C:\mitsuba\>rename dependencies_windows dependencies</code>

Noted that some codes in this tool are from Mitsuba source code and Mitsuba is licensed under the terms of Version 3 of the GNU General Public License, so if you use this code in your project, please follow the rules described in the license.
