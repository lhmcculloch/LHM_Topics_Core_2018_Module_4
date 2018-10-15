#LHM
#ATiB Module 4 Code
#10/14/18

#Objective 1: Analyze overlap of K562 and GM12878 GROcap data from Core et al. with ENCODE DNase hypersensitive sites
#Note: A bunch of this was done in as command line code in Terminal, but I'm including it in the same file as my R code for ease of access.
#Commenting it out, but including it in this file.

#Obtained GROcap TSS data from supplementary materials: tss_all_gm12878.bed and tss_all_k562.bed
#Obtained DHS data from ENCODE:
#Non-cell-line specific: wgEncodeRegDnaseClusteredV3.bed.gz
#Cell-line-specific: Used narrowPeak and broadPeak replicate 1 for each cell line

'

#Sorting the TSS data:
sort-bed tss_all_k562.bed > tss_all_k562_sorted.bed
sort-bed tss_all_gm12878.bed > tss_all_gm12878_sorted.bed


#Sorting the DHS data:
sort-bed wgEncodeRegDnaseClusteredV3.bed > wg_DHS_hg19_sorted.bed
sort-bed K562_narrowPeak_rep1_hg19.bed > K562_NP_hg19_sorted.bed
sort-bed GM12878_narrowPeaks_rep1_hg19.bed > GM12878_NP_hg19_sorted.bed
sort-bed K562_broadPeaks_hg19.bed > K562_BP_hg19_sorted.bed
sort-bed GM12878_broadPeaks_hg19.bed > GM12878_BP_hg19_sorted.bed


#Overall TSS counts:
wc -l tss_all_k562_sorted.bed #128471 --> use this as divisor for K562 results
wc -l tss_all_gm12878_sorted.bed #117613 --> use this as divisor for GM12878 results


#DHS starting counts:
wc -l wg_DHS_hg19_sorted.bed #Full DHS list #1867665
wc -l GM12878_NP_hg19_sorted.bed #narrowPeak DHS list #56393
wc -l GM12878_BP_hg19_sorted.bed #broadPeak DHS list #51916
wc -l K562_NP_hg19_sorted.bed #narrowPeak DHS list #169300
wc -l K562_BP_hg19_sorted.bed #broadPeak DHS list #156253

#Using bedops --element-of to see which TSSs align to at least one DHS. Looks at each TSSs (1 bp each) and sees
#if each one overlaps a DHS. Calculate (TSS reads overlapping a DHS) / (total TSS reads) to get fraction overlap.

#Using the full DHS list
bedops --element-of 1 tss_all_k562_sorted.bed wg_DHS_hg19_sorted.bed > K562_to_DHS_all_hg19.bed
bedops --element-of 1 tss_all_gm12878_sorted.bed wg_DHS_hg19_sorted.bed > GM12878_to_DHS_all_hg19.bed
wc -l K562_to_DHS_all_hg19.bed #118494 --> 92.2%
wc -l GM12878_to_DHS_all_hg19.bed #112436 --> 95.6%

#narrowPeak cell-line-specific data
bedops --element-of 1 tss_all_k562_sorted.bed K562_NP_hg19_sorted.bed > K562_to_DHS_NP_hg19.bed
bedops --element-of 1 tss_all_gm12878_sorted.bed GM12878_NP_hg19_sorted.bed > GM12878_to_DHS_NP_hg19.bed
wc -l K562_to_DHS_NP_hg19.bed #65329 --> 50.9%
wc -l GM12878_to_DHS_NP_hg19.bed #37299 --> 31.7%

#broadPeak cell-line-specific data
bedops --element-of 1 tss_all_k562_sorted.bed K562_BP_hg19_sorted.bed > K562_to_DHS_BP_hg19.bed
bedops --element-of 1 tss_all_gm12878_sorted.bed GM12878_BP_hg19_sorted.bed > GM12878_to_DHS_BP_hg19.bed
wc -l K562_to_DHS_BP_hg19.bed #106972 --> 83.3%
wc -l GM12878_to_DHS_BP_hg19.bed #62924 --> 53.5%

#Looking at intersect results in bedops for a more global picture of how TSSs and DHSs overlap:

#GM12878 intersection results (for Venn diagrams):
#All DHSs:
bedops --intersect tss_all_gm12878_sorted.bed wg_DHS_hg19_sorted.bed > GM12878_to_DHS_all_hg19_intersect.bed
wc -l GM12878_to_DHS_all_hg19_intersect.bed #102538

#NP DHSs:
bedops --intersect tss_all_gm12878_sorted.bed GM12878_NP_hg19_sorted.bed > GM12878_to_DHS_NP_hg19_intersect.bed
wc -l GM12878_to_DHS_NP_hg19_intersect.bed #33068

#BP DHSs:
bedops --intersect tss_all_gm12878_sorted.bed GM12878_BP_hg19_sorted.bed > GM12878_to_DHS_BP_hg19_intersect.bed
wc -l GM12878_to_DHS_BP_hg19_intersect.bed #56240


#K562 intersection results (for Venn diagrams):
#All DHSs:
bedops --intersect tss_all_K562_sorted.bed wg_DHS_hg19_sorted.bed > K562_to_DHS_all_hg19_intersect.bed
wc -l K562_to_DHS_all_hg19_intersect.bed #107766

#NP DHSs:
bedops --intersect tss_all_K562_sorted.bed K562_NP_hg19_sorted.bed > K562_to_DHS_NP_hg19_intersect.bed
wc -l K562_to_DHS_NP_hg19_intersect.bed #58237

#BP DHSs:
bedops --intersect tss_all_K562_sorted.bed K562_BP_hg19_sorted.bed > K562_to_DHS_BP_hg19_intersect.bed
wc -l K562_to_DHS_BP_hg19_intersect.bed #95785
'



#Also did some trial work on figure 2a. Several of us worked on this in parallel, and we didn't use the code below
#in our final version, but I'm including it for reference. We also realized partway along that we needed the
#protein-coding gene data, not the transcript data; this version is from before that conclusion, and sorts for the
#transcript data. For the final version that generated the figure, see Runyu's code.



#Started with the Gencode annotated gene/transcript list, called gencode.bed
#Pulling out the portions we initially thought we were interested in. Undecided about whether to use +, -, or both at this point
#awk '($8=="transcript")' gencode.bed > gentxn.bed
#awk '($6=="+")' gentxn.bed > gentxn_plus.bed
#awk '($6=="-")' gentxn.bed > gentxn_minus.bed

#Pulling out transcripts > 9 kb long, so we could take 3 x 3 kb segments without overlaps
#awk '($3 - $2 >= 9000)' gentxn_plus.bed > gtplus9000.bed #51023 lines remain at this point, of 99306 originally
#awk '($3 - $2 >= 9000)' gentxn_minus.bed > gtminus9000.bed #49564 of 97214 remain
#awk '($3 - $2 >= 9000)' gentxn.bed > gt9000.bed

#Moved the files into R at this point for further manipulation:
gt9000_bed <- read.table("gt9000.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
gtplus9000_bed <- read.table("gtplus9000.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
gtminus9000_bed <- read.table("gtminus9000.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

library(dplyr)

#Getting listed coordinates for first 3 kb, middle 3 kb, and last 3 kb:

#First 3 kb:
tss3000_all <- select(gt9000_bed, V1, V2, V3) #Pulling out the columns that I'm interested in. Starting with this for the combined + and - list, but could do these separately.
for (i in 1:length(rownames(tss3000_all))) {
  tss3000_all$V3[i] <- (tss3000_all$V2[i] + 3000)
}

#Middle 3 kb:
mid3000_all <- select(gt9000_bed, V1, V2, V3)
for (i in 1:length(rownames(end3000_all))) {
  midpt = as.integer((mid3000_all$V2[i] + mid3000_all$V3[i])/2)
  mid3000_all$V2[i] <- (midpt - 1500)
  mid3000_all$V3[i] <- (midpt + 1500)
}

#Last 3 kb:
end3000_all <- select(gt9000_bed, V1, V2, V3)
for (i in 1:length(rownames(end3000_all))) {
  end3000_all$V2[i] <- (end3000_all$V3[i]-3000)
}

#Next, wanted to figure out where the GRO-cap reads overlapped each range. For each overlap, wanted to track a
#count of how many alignments occurred at each position (1-3000), and keep a score for each position (for potential
#use later in normalization). Opted to try to keep both an absolute value score and a score using +/- values for +/- strand, respectively.
#Ultimately tried several things for this approach, including a couple of loops (too slow) and bedmap. Listing for
#reference.

#Adding in the bed files:
GM12878_gro_all_bed <- read.table("GM12878_gro_all.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
K562_gro_all_bed <- read.table("K562_gro_all.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

#Setting up scoring
tss3000_score_df <- data.frame("Position" = 1:3000, "Count" = 0, "Abs_Score" = 0, "Signed_Score" = 0)

#Tried a few variants of this loop with some breaks, additional conditional statements, etc., to try to speed things
#up, but ultimately concluded that it wasn't feasible and abandoned it.
for (i in 1:length(rownames(GM12878_gro_all_bed))) {

  for (j in 1:length(rownames(tss3000_all))) {
 
    if ((GM12878_gro_all_bed$V1[i]==tss3000_all$V1[j]) && (GM12878_gro_all_bed$V2[i] >= tss3000_all$V2[j]) && (GM12878_gro_all_bed$V2[i] < tss3000_all$V3[j])) {

      for (k in 1:3000) {
        if ((tss3000_all$V2[j] + k - 1) == GM12878_gro_all_bed$V2[i]) {#If a position is a match
          tss3000_score_df$Count[k] <- (tss3000_score_df$Count[k] + 1) #Add one to the count
          tss3000_score_df$Abs_Score[k] <- (tss3000_score_df$Abs_Score[k] + abs(GM12878_gro_all_bed$V5[i]))
          tss3000_score_df$Signed_Score[k] <- (tss3000_score_df$Signed_Score[k] + GM12878_gro_all_bed$V5[i])

        }
      }
    }
  }
}

#An alternative loop. Again, too slow, so abandoned it:
tss3000_score_alt_df <- data.frame("Position" = 1:3000, "Count" = 0, "Abs_Score" = 0, "Signed_Score" = 0)

for (i in 1:length(rownames(GM12878_gro_all_bed))) {
  for (j in 1:length(rownames(tss3000_all))) {
    if ((GM12878_gro_all_bed$V1[i]==tss3000_all$V1[j]) && (GM12878_gro_all_bed$V2[i] >= tss3000_all$V2[j]) && (GM12878_gro_all_bed$V2[i] < tss3000_all$V3[j])) { #chromosomes match and the coordinate is between the two other coordinates, i.e., it's a match
      k <- (tss3000_all$V2[j] - GM12878_gro_all_bed$V2[i] + 1)
      tss3000_score_alt_df$Count[k] <- (tss3000_score_alt_df$Count[k] + 1) #Add one to the count
      tss3000_score_alt_df$Abs_Score[k] <- (tss3000_score_alt_df$Abs_Score[k] + abs(GM12878_gro_all_bed$V5[i]))
      tss3000_score_alt_df$Signed_Score[k] <- (tss3000_score_alt_df$Signed_Score[k] + GM12878_gro_all_bed$V5[i])
    }
  }
}


#Contrary to what the documentation suggested, another group member said that bedmap map-range could give the entirety
#of each portion of the overlapping sequence, not just the overlap coordinates. #Try using bedmap map-range. It doesn't
#just give the overlap ranges, it gives the entire thing. If I run it with my samples and my 3000 bp regions, I
#think it should give me the info I want. Then I can use subtraction to get the relative position information.

#For this approach, decided just to use the + strand data.

tss3000_plus <- select(gtplus9000_bed, V1, V2, V3) #Pulling out the columns that I'm interested in.

#Pulling our regions of interest
for (i in 1:length(rownames(tss3000_plus))) {
  tss3000_plus$V3[i] <- (tss3000_plus$V2[i] + 3000)
}

write.table(tss3000_plus, file = "~/Desktop/Topics_M4/tss3000_plus.bed", row.names = FALSE, col.names = FALSE) #Got the file out

#Back to Terminal:

#sort-bed tss3000_plus.bed > tss3000_plus_sorted.bed
#sort-bed GM12878_gro_plus.bed > GM_gro_plus_sorted.bed
#bedmap --delim '\t' --echo-map-range --echo GM_gro_plus_sorted.bed tss3000_plus_sorted.bed > GM_to_TSS_3000_overlaps.txt
#Didn't produce expected overlaps, but when another group member did the same thing, it worked.
#Group wasn't clear what the issue was, but proceeded from here with the working overlap data.

