<br/><br/>

<br/><br/>

<br/><br/>

<p style="font-family: &#39;Fira Mono&#39;, monospace; font-size: 2.5rem;">
**myTAI** unlocks statistically-informed analysis of evolutionary
signals hidden in the gene expression data (transcriptome)
</p>

<br/><br/>

<br/><br/>

<br/><br/>

    library(myTAI)
    # obtain an example phylo-expression object
    data("example_phyex_set")
    # plot away!
    myTAI::plot_signature(example_phyex_set)  

<img src="man/figures/plot_signature.svg" alt="plot_signature" style="display: block; margin: auto;" />

<br/><br/>

<br/><br/>

<br/><br/>

<br/><br/>

<br/><br/>

<p style="font-family: &#39;Fira Mono&#39;, monospace; font-size: 1.5rem;">
Initially inspired by the need to statistically evaluate the ‘molecular
hourglass’ pattern that linked transcriptome evolution to 19th-century
embryo studies, the original myTAI package has since evolved to
facilitate fast and robust analysis of *any* transcriptome evolution
patterns, from embryo development[1] to cancer progression and drug
perturbation experiments.
</p>

<br/><br/>

<br/><br/>

------------------------------------------------------------------------

## What can `myTAIv2` do?

<p style="font-family: &#39;Fira Mono&#39;, monospace; font-size: 1.5rem;">
⚡ `Ultra-fast computations` → rapid permutation tests that scales to
single-cell genomics.  
🦾 `Moved on from traditional R data structures to modern S7 classes` →
optimised handling of large ‘phylo-expression’ datasets.  
🔨 `Break your hourglass` → uncover the genes driving your transcriptome
evolution patterns.  
⌛ `Aesthetic plots` → publication-ready information-rich figures!
</p>

<br/><br/>

<p style="font-family: &#39;Fira Mono&#39;, monospace; font-size: 2rem;">
Get started with myTAI →
<a href="articles/myTAI.html" class="btn btn-outline-light" style="font-family: 'Fira Mono', monospace; font-size: 2rem;">
🍹 </a>
</p>

------------------------------------------------------------------------

<br/><br/>

<br/><br/>

<br/><br/>

<br/><br/>

<p style="font-family: &#39;Fira Mono&#39;, monospace; font-size: 1.5rem;">
Cite myTAI for your own research.
</p>
<p style="font-family: &#39;Fira Mono&#39;, monospace; font-size: 1.5rem;">

> Drost et al. **myTAI: evolutionary transcriptomics with R**.
> *Bioinformatics* 2018, 34 (9), 1589-1590.
> [doi:10.1093](https://doi.org/10.1093/bioinformatics/btx835)

</p>
<p style="font-family: &#39;Fira Mono&#39;, monospace; font-size: 1.5rem;">
And follow [`@DrostLab`](https://drostlab.com/) for more bioinformatics
and digital biology software solutions.
</p>

<br/><br/>

<br/><br/>

[1] some early studies include [Kalinka et al. (2010)
Nature](https://www.nature.com/articles/nature09634), [Domazet-Lošo &
Tautz (2010) Nature](https://www.nature.com/articles/nature09632),
[Quint et al. (2012)
Nature](https://www.nature.com/articles/nature11394)
