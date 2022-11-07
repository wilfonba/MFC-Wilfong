<p align="center">
  <img src="doc/res/banner.png" alt="MFC Banner" width="500"/></center>
</p>

<p align="center">
  <a href="http://dx.doi.org/10.1016/j.cpc.2020.107396">
    <img src="https://zenodo.org/badge/doi/10.1016/j.cpc.2020.107396.svg" />
  </a>
  <a href="https://github.com/MFlowCode/MFC/actions">
    <img src="https://github.com/MFlowCode/MFC/actions/workflows/ci.yml/badge.svg" />
  </a>
  <a href="https://github.com/MFlowCode/MFC/actions">
    <img src="https://github.com/MFlowCode/MFC/actions/workflows/doc.yml/badge.svg" />
  </a>
  <a href="https://lbesson.mit-license.org/">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg" />
  </a>
  <a href="https://github.com/MFlowCode/MFC/commit/">
    <img src="https://badgen.net/github/last-commit/MFlowCode/MFC" />
  </a>
  <a href="https://hub.docker.com/repository/docker/henryleberre/mfc">
    <img src="https://shields.io/docker/image-size/henryleberre/mfc" />
  </a>
  <a href="https://zenodo.org/badge/latestdoi/198475661">
    <img src="https://zenodo.org/badge/198475661.svg" alt="DOI">
  </a>
</p>

Welcome to the home of MFC!
MFC simulates compressible multi-component and multi-phase flows, amongst other things. 
It scales ideally to tens of thousands of GPUs.
MFC is a Fortran codebase that makes use of metaprogramming to keep the codebase readable and compact.
Please contact the developers, like [Spencer](mailto:shb@gatech.edu), if you have any questions.
We have an active Slack channel to help ease in new MFC users and support development.

MFC has both API and high-level documentation. 
It is available on [the website](https://mflowcode.github.io/) and in markdown format at [doc/landing/readme.md](doc/landing/readme.md).

If you use MFC, consider citing it:
* <a href="https://doi.org/10.1016/j.cpc.2020.107396">
    S. H. Bryngelson, K. Schmidmayer, V. Coralic, K. Maeda, J. Meng, T. Colonius (2021) Computer Physics Communications 4655, 107396
</a>

```
@article{Bryngelson_2021,
  title = {MFC: An open-source high-order multi-component, multi-phase, and multi-scale compressible flow solver},
  author = {Spencer H. Bryngelson and Kevin Schmidmayer and Vedran Coralic and Jomela C. Meng and Kazuki Maeda and Tim Colonius},
  journal = {Computer Physics Communications},
  doi = {10.1016/j.cpc.2020.107396},
  year = 2021,
  pages = {107396},
}
```

## License
 
Copyright 2022.
MFC is under the MIT license (see [LICENSE](LICENSE) file for full text).

## Acknowledgements
 
<p align="justify">
   MFC development was supported by multiple current and past grants from the US Office of Naval Research (ONR), the US National Institute of Health (NIH), and the US National Science Foundation (NSF).
  MFC computations utilize the Extreme Science and Engineering Discovery Environment (XSEDE, now ACCESS), under allocations TG-CTS120005 (PI Colonius) and TG-PHY210084 (PI Bryngelson) and OLCF Summit under allocation CFD154 (PI Bryngelson).
</p>
