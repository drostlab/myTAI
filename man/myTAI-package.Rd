\name{myTAI-package}
\alias{myTAI-package}
\alias{myTAI}
\docType{package}
\title{
A package to perform phylotranscriptomics analyses for Evolutionary Developmental Biology research.
}
\description{

\emph{Phylotranscriptomics} defines the concept of combining genetic sequence information with gene expression levels to capture evolutionary signals during developmental processes.

Common measures in \emph{phylotranscriptomics} analyses are:

\emph{The Transcriptome Age Index:}

The phylogenetically based transcriptome age index (TAI) measures the mean transcriptome age for a given developmental time course (Domazet-Loso & Tautz, 2010). 


\emph{The Transcriptome Divergence Index:}

The sequence divergence based transcriptome divergence index (TDI) quantifies the mean transcriptome divergence between two different species for a given developmental time course (Quint et al., 2012).


This package provides a broad toolbox for common phylotranscriptomics analyses and also
includes statistical tests to verify observed phenomena and transcriptomic patterns.

Using the \emph{Rcpp} package, all computationally expansive functions have been written in C++ to enable real-time analytics on phylotranscriptomics datasets.

In the \emph{Examples} section a common phylotranscriptomics workflow using the \emph{myTAI} package is presented and can be adapted to any standardized phylotranscriptomics dataset.

\bold{\emph{A short history of phylotranscriptomics:}}

\emph{Phylotranscriptomics} arose from the interdisciplinary field today known as \emph{Evolutionary Developmental Biology} which combines research in developmental biology with evolutionary biology to study the evolution of developmental processes and the evolvability and diversity of complex life on earth. A great review and introduction to \emph{Evolutionary Developmental Biology} can be found in (Sean B. Carroll, 2008).

In his famous book \emph{On the Origin of Species} Charles Darwin writes: 

\emph{It has been disputed at what period of life the causes of variability, whatever they may be, generally act; whether during the early or late period of development of the embryo, or at the instant of conception.} (Charles Darwin, 1859)[p. 9]. 

Karl Ernst von Baer was the first person linked with the observation that animal embryos morphologically resemble each other during a specific period of embryo development (Karl Ernst von Baer, 1828). 

Astonished by his observations, von Baer wrote a letter to Charles Darwin:

\emph{In my possession are two little embryos in spirit [alcohol], whose names I have omitted to attach, and at present I am quite unable to say to what class they bolong. They may be lizards or small birds, or very young mammalia, so complete is the similarity in the mode of formation of the head and trunk in these animals.} (Jerry A. Coyne, 2009)[p. 78].

Darwin later writes in his On the Origin of Species:

\emph{The points of structure, in which the embryos of widely different animals of the same class resemble each other, often have no direct relation to their conditions of existence. We cannot, for instance, suppose that in the embryos of the vertebrata the peculiar loop-like course of the arteries near the branchial slits are related to familiar conditions - in the young mammal which is nourished in the womb of its mother, in the egg of the bird which is hatched in a nest, and in the spawn of a frog under water.} (Jerry A. Coyne, 2009)[p. 82], (Charles Darwin, 1859).


These quotes from \emph{von Baer} and \emph{Charles Darwin} shall illustrate the long fascination of biologists for this mysterious period during animal embryogenesis in which animals of different phyla morphologically resemble each other observably.

Different hypotheses have been proposed to explain this mysterious morphological series of embryo development. The most commonly cited models are: 

(1) the \emph{Early Conservation Model} which proposes that the early developmental stages are foundational and that morphological complexity is increasing during embryogenesis; 

(2) the \emph{Developmental Hourglass Model} reflecting the evolutionary most conserved stage during embryogenesis and is believed to be the stage at which coordination between growth and patterning is regulated to build the body-plan; 

(3) the \emph{Adaptive Penetrance Model} proposing that mid-embryogenesis stages are the most important stages for beneficial mutations; and 

(4) the \emph{Ontogenetic Adjacency Model} proposing that small changes are most likely between events that are adjacent in the developmental sequence of events (Alex T. Kalinka and Pavel Tomancak, 2012).


Different studies provide convincing evidence for a developmental hourglass model, but morphological and molecular phenomena only predicted by the \emph{early conservation model} causes an ongoing depate about the developmental and evolutionary processes that are shaping this morphological transition (Katherine E. Willmore, 2012).

The field today referred to as \emph{phylotranscriptomics} arose from this debate and provides a new perspective on this debate by combining techniques that quantify the evolutionary relationships on the molecular level with techniques that quantify the molecular control mechanisms of development.


K. Sander introduced the term \emph{phylotypic stage} to account for the striking morphological similarity of embryos from different phylogenetic groups during this specific period of embryogenesis (K. Sander, 1983). The term \emph{developmental hourglass} was used by Denis Duboule and Rudolf A. Raff to illustrate the morphological series of embryo resemblance during embryogenesis (Denis Duboule, 1994),  (Rudolf A. Raff, 1996).


Striking results from recent studies propose that the developmental hourglass model can be translated to the molecular level, and furthermore, also follow a molecular hourglass pattern (Quint et al., 2012). 

O. C. Ivanov was the first person who tested von Baer’s and Haeckel’s laws of embryo morphology on the molecular level based on mammalian globin sequence data (Hazkani-Covo et al., 2005) and thus provided a first idea how to translate the morphological phenomena into molecular biology.

Ivanov concluded that Haeckel’s biogenetic law also exists on the molecular level based on his observations concerning the sequence identities of embryonic hemoglobin chains (O. C. Ivanov, 1987). Duboule who proposed the developmental hourglass model suggests that during the phylotypic stage vertebrate Hox genes have a specific gene expression differing from the gene expression during early and late embryogenesis.

Duboule additionally proposed that vertebrate developmental strategies converge towards the phylotypic stage to establish a stable bodyplan. For Duboule the phylotypic stage reflects the stage or period in which Hox gene constraints are most stringent. Hox genes are known to be involved in axis-formation during vertebrate development and therefore, Duboule concluded that the phylotypic stage or phylotypic period can be seen as a stage of reduced evolutionary allowance for modifications in axis-formation (Denis Duboule, 1994).

Furthermore, Raff concluded that during the phylotypic period major gene regulators governing fundamental processes of organ formation are activated, whereas no organ formation has to be governed during early embryo development. During late embryo development only very specialized organs are formed which are governed by intermediate gene regulators (Rudolf A. Raff, 1996).

Later, Hazkani-Covo at al. found that more evolutionarily conserved genes are expressed during the phylotypic stage of mouse embryo development than it is the case during early and late embryo development. Based on gene expression studies of orthologous genes in human and mouse, they assume that orthologous genes expressed in the phylotypic stage should resemble each other more closely than orthologs expressed in early and late embryogenesis.

Additionally, Hazkani-Covo at al. assume that a particular gene that is expressed during a specific developmental stage in the mouse is also expressed at the same stage in human. Most of their findings were statistically not significant and therefore, they concluded that there is some evidence for the existence of a molecular hourglass model during mouse embryogenesis, but a true confirmation can only be done when more data is available (Hazkani-Covo et al., 2005). 

Five years later, Kalinka et al. used a similar approach to test for the existence of a molecular hourglass pattern based on gene expression data from six Drosophila species with an evolutionary distance up to 40 Ma (Alex T. Kalinka et al., 2010). 

They found that gene expression of the 3,019 orthologous genes known to be expressed during Drosophila embryogenesis is most conserved during the arthropod phylotypic stage and that selective constraints reach their maximum during the phylotypic stage. They also found that gene expression is more resistant to evolutionary changes during the phylotypic stage than during early and late embryogenesis. Kalinka et al. conclude that in accordance with Raff, an increase in global genetic interactions constrains evolutionary changes within the phylotypic period (Alex T. Kalinka et al., 2010).

In the same year Domazet-Loso and Tautz proposed the existence of transcriptomic hourglass patterns in zebrafish, fly, mosquito, and worm. Based on phylotranscriptomics studies, Domazet-Loso and Tautz argue that \emph{the transcriptomic hourglass pattern highlights a link between evolutionary innovations and the emergence of novel genes.} As a result, they suggest a new definition of the phylotypic stage as \emph{[...] the ontogenetic progression during which the oldest gene set is expressed, either because this is the phase with the lowest opportunity for lineage-specific adaptations, or because it is internally so constrained that newly evolved genes cannot become integrated.} (Tomislav Domazet-Loso and Diethard Tautz, 2010).


Today, the molecular hourglass pattern has also been observed in nematodes, frogs, and turtles (Alex T. Kalinka and Pavel Tomancak, 2012). Piasecka et al. used a modularization approach by decomposing genes into independent stage-specific modules to study the Early Conservation Model and the Developmental Hourglass Model in zebrafish. Piasecka et al. found that gene duplication and de novo gene emergence is under strong constraints during early development, whereas regulatory elements have been found to be most conserved for TFs expressed at mid-embryogenesis. However, they did not observe an hourglass pattern for gene sequences, gene age, and gene expression using their modularization approach (Piasecka et al., 2013). These findings mark the ongoing debate about the existence of an \emph{Early Conservation Model} or \emph{Hourglass Model} also on the molecular level. 

Alicia N. Schep and Boris Adryan also found that TF expression mostly fits with the hourglass model of development. They show for vertebrates and arthropods, that animal-specific TF expression is most similar between two species during mid-embryogenesis and conclude that their results support the hypothesis raised by (Kalinka et al., 2010) that although species pass through similar developmental events, the regulatory re-wiring by these events can diverge extensively (Alicia N. Schep and Boris Adryan, 2013).

Despite the debate, whether a developmental hourglass exists on the morphological level, the molecular level, or on neither of both levels (Alex T. Kalinka and Pavel Tomancak, 2012), studying the connections between developmental processes, their genetic regulations, and their morphological outcome, have shown to provide fascinating facts of life strategies and phylogenetic relationships. Most of the studies show that animal-specific TFs and cis-regulatory elements share more evolutionary similarity during mid-embryogenesis than during early and late development.


}
\details{
\tabular{ll}{
Package: \tab myTAI\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-7\cr
License: \tab >= GPL-2\cr
}

}
\author{
Hajk-Georg Drost

Maintainer: Hajk-Georg Drost <hajk-georg.drost@informatik.uni-halle.de>}
\references{

Baer von KE. 1828. "Entwicklungsgeschichte der Thiere: Beobachtung und Re-
flexion". (Vol. 1). Koenigsberg, Borntraeger.

Carroll SB. 2008. "Evo-devo and an expanding evolutionary synthesis: A genetic theory of morphological evolution". Cell (134): 25–36.

Coyne JA. 2009. "Why Evolution Is True". Oxford University Press.

Darwin C. 2008. "On the Origin of Species". Oxford World’s Classics, Oxford University Press, Great Clarendon Street, Oxford.

Domazet-Loso T and Tautz D. 2010. "A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns". Nature (468): 815-818.

Drost HG. 2013. "A bioinformatics approach to study the origin of embryogenesis in plants and animals". Master's Thesis, Martin-Luther University University Halle-Wittenberg.

Duboule D. 1994. "Temporal colinearity and the phylotypic progression: a basis for the stability of avertebrate bauplan and the evolution of morphologies through heterochrony". Development, pages 135–142.

Eddelbuettel D and Francois R. 2011 "Rcpp: Seamless R and C++ Integration". Journal of
Statistical Software, 40(8), 1-18. URL http://www.jstatsoft.org/v40/i08/.

Eddelbuettel D. 2013 "Seamless R and C++ Integration with Rcpp". Springer, New York. ISBN 978-1-4614-6867-7.

Hazkani-Covo E et al. 2005. "In search of the vertebrate phylotypic stage: A molecular examination of the developmental hourglass model and von baer’s third law". JOURNAL OF EXPERIMENTAL ZOOLOGY (304): 150– 158.

Ivanov OC. 1987. "New evidence for validity of haeckel’s law on molecular level". Naturwissenschaften (74): 40–42.

Kalinka AT and Tomancak P. 2012. "The evolution of early animal embryos: conservation or divergence?". Trends in Ecology and Evolution (27): 385–393.

Kalinka A et al. 2010. "Gene expression divergence recapitulates the develop- mental hourglass model". Nature (468) :811–814.

Katherine KE. 2012. "The Body Plan Concept and Its Centrality in Evo-Devo". Evo Edu Outreach (5): 219–230.

Piasecka B et al. "The hourglass and the early conservation models - co-existing evolutionary patterns in vertebrate development". PLOS Genetics 9(4): e1003476.

Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.

R Core Team. 2013. "R: A language and environment for statistical computing. R Foundation for Statistical Computing". R Foundation for Statistical Computing, Vienna, Austria.

Schep AN and Adryan B. 2013. "A Comparative Analysis of Transcription Factor Expression during Metazoan Embryonic Development".  PLOS ONE (8): e66826.


}
\keyword{ phylotranscriptomics }
\seealso{

}

