# What do we see when we look at networks?
This is the repository for the data, images and code of the paper "What do we see when we look at networks" by Tommaso Venturini, Mathieu Jacomy, Pablo Jensen".

* networks: contains the GEXF files of the networks presented in the paper.
* HiRes images: contains the hi-resolution version of all the images created for the paper
* Supplementary_Clustering_Comparison.pdf: extend the comparison between different clustering algorithms discussed in section 4
* Network images: those images were produced by applying a script to the different networks. Note: for Karate_FA_Default and Karate_FA_LinLog_grav0 the images produced have different settings (bigger nodes) to improve the paper's readability.
* Graph recipes scripts: script used to compute and visualize different statistics on a network.

---
### How to use a graph recipes script

1. Go on the following URL: http://tools.medialab.sciences-po.fr/graph-recipes/
2. Upload the GEXF network of your choice
3. Pick an empty recipe ("Nothing, I'm fine")
4. Open the script in a text editor and copy-paste its content in the Javascrip area
5. Run the recipe to produce the computations and images
