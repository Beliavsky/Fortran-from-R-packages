
---
title: "Credit Default Swaps with R"
tags:
- credit default swap
- pricing
- R package
date: "7 January 2025"
output: pdf_document
authors:
- name: Yanrong Song
  orcid: "0009-0003-5738-8384"
  email: yrsong129@gmail.com
  affiliation: 1
- name: Zijie Zhu
  email: zijie.miller.zhu@gmail.com
  affiliation: 1
- name: David Kane
  email: david.kane@gmail.com
  affiliation: 2
- name: Heidi Chen
  email: s.heidi.chen@gmail.com
  affiliation: 2
- name: Yuanchu Dang
  email: yuanchu.dang@gmail.com
  affiliation: 1
- name: Yang Lu
  email: yang.lu2014@gmail.com
  affiliation: 3
- name: Kanishka Malik
  email: kanishkamalik@gmail.com
  affiliation: 3
- name: Skylar Smith
  email: sws2@williams.edu
  affiliation: 3
bibliography: paper.bib
affiliations:
- name: Columbia University, United States
  index: 1
- name: Harvard University, United States
  index: 2
- name: Williams College, United States
  index: 3
---



# Summary

The `creditr` \textbf{R} package is an open-source tool for modeling and valuing credit default swaps (CDS). CDS are financial derivatives that allow investors to hedge or speculate on credit risk by providing protection against the default of a third party [@creditDerivatives]. Widely used in financial markets, CDS require complex pricing models to account for variables such as recovery rate, coupon, interest rate, and default probabilities. The `creditr` package  provides users with tools to calculate CDS spreads, upfront payments, and other related quantities using the *ISDA Standard Model* [@isdaCDS2022], the industry benchmark for CDS pricing. The package enables users to replicate outputs from proprietary platforms such as Bloomberg and Markit [@rates], offering a transparent and scriptable alternative for CDS pricing and analytics. Designed for fixed income researchers and finance practitioners, `creditr` makes quantitative credit modeling accessible within the \textbf{R} programming ecosystem.


# Statement of need

Despite the importance of credit derivatives in modern financial markets, most CDS pricing tools remain embedded in proprietary platforms such as Bloomberg or Markit and have limited accessibility and reproducibility outside of these platforms, imposing significant entry barrier to new practitioners and academics. Existing open-source financial libraries also have very limited support for credit derivatives, especially those following the ISDA Standard Model.

The `creditr` package addresses this gap by providing a fully open-source implementation of key CDS pricing routines. Suppose a researcher wishes to compute the fair upfront cost of buying CDS protection on a corporate bond with a 5-year maturity and a quoted spread of 150 basis points. Using `creditr`, the user can specify the notional amount, recovery rate, interest rate term structure, and default probability curve, then call the `calc_upfront()` function to obtain the fair premium or cost. The results can be compared directly with Bloomberg's CDS Calculator output, ensuring practical consistency with industry tools.

Overall, `creditr` enables users to:

- Compute the fair spread or upfront payment for a CDS given market inputs.
- Perform sensitivity analysis with respect to interest rates, recovery rates, and default curves.
- Replicate the outputs of the Bloomberg CDS Calculator and Markit CDS pricers and demystify the mechanics behind CDS pricing.

Please note that the **creditr** \textbf{R} package is a wrapper around ISDA's standard CDS model in \textbf{C++}.

# References