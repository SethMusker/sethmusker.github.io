---
title: "test formulae"
date: 2025-11-26
categories: []
tags: []
output:
  md_document:
    variant: markdown
    preserve_yaml: true
math: true
---

# test formulae

## Introduction

Here is the model equation:

$$
\mathrm{rate} =
\frac{
  r_{\mathrm{tref}}\;\exp\!\Bigl(-\frac{e}{k}\Bigl(\frac{1}{T + 273.15} -
  \frac{1}{T_{\mathrm{ref}} + 273.15}\Bigr)\Bigr)
}{
  1 + \left(\frac{e}{e_h - e}\right)\;
  \exp\!\Bigl(\frac{e_h}{k}\Bigl(\frac{1}{T_{\mathrm{opt}} + 273.15} -
  \frac{1}{T + 273.15}\Bigr)\Bigr)
}
$$
