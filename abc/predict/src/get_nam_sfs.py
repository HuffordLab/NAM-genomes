from nam_sfs_functions import *
import argparse

parser = argparse.ArgumentParser(description="Generate SFS for nam founders according to regions provided in a bedfile.")

parser.add_argument("-b", "--bedfile", type=str, required = True,
            help = "Name/path of bedfile giving intervals for sfs.")

parser.add_argument("-v", "--variantfile", type=str, required = True,
            help = "Name/path of variant file --either vcf or a structural variant file.")


parser.add_argument("-t", "--variant_type", type=str, required = True, choices = {"vcf", "sv"},
            help = "variant input type based on what type variantfile (vcf or sv) is provided.")

parser.add_argument("-T", "--sv_type", type=str, default = None,
            help = "When structural variant file is input, this filters to only comput sfs for chosen kind. If absent, all are computed. Multiple can be provided, assumes name is broken by and underscore, without spaces")


#parser.add_argument('-a', action="store_true", default=False)

parser.add_argument("-A", "--any_site", action='store_false', default = True,
                help = "Switch allow any sites in the VCF to be called, rather than just biallelic ones.")

parser.add_argument("-f", "--foldfile", type=str,
            help = "Optional Name/path of file giving sites that are 0 or 4 fold. input lines should be chrom\t1 indxed position\t[04]")

parser.add_argument("-F", "--foldtype", type=str,
            help = "If foldfile is provided, gives 0 or 4, for selecting 0 or 4 fold sites from those listed.")

parser.add_argument("-I", "--sfs_size", type = int,
            help = "Total number of categories for the SFS excluding invariant sites category. Folded, so the end of the array will be zeros. Default is 26.")

parser.add_argument("-H", "--heterozygous_prop", type = float, default = 0,
                    help = "Proportion of individuals allowed to be heterozygous. Defaults to 0.")

parser.add_argument("-M", "--missing_prop", type = float, default = 0,
                    help = "Proportion of individuals allowed to be missing. Defaults to 0.")

args = parser.parse_args()


bed_df = pd.read_csv(args.bedfile, sep = "\t", names = ['chrome','start','end'])


if args.sv_type:
    sv_list = args.sv_type.split("_")
else:
    sv_list = None

if args.foldfile:
    fold_dict = make_site_dict(args.foldfile, args.foldtype)
    SFS = build_sfsdf(bed_df, sfs_size = args.sfs_size)
    SFS = file_parse(args.variantfile, args.variant_type, bed_df, SFS,
                     i_count = args.sfs_size, 
                     missing_prop = args.missing_prop, bi_allelic = args.any_site, 
                     heterozygous_prop = args.heterozygous_prop, 
                     site_dict = fold_dict, sv_type = sv_list
                    )
else:
    SFS = build_sfsdf(bed_df, sfs_size = args.sfs_size)
    SFS = file_parse(args.variantfile, args.variant_type, bed_df, SFS,  
                     i_count = args.sfs_size, bi_allelic = args.any_site,
                     missing_prop = args.missing_prop, 
                     heterozygous_prop = args.heterozygous_prop, 
                     sv_type = sv_list
                    )

for idx, sfs in SFS.items():
    bed_row = bed_df.iloc[idx, :]
    print('\t'.join([str(i) for i in sfs]))

#bed_df = pd.read_csv("data/ref/Zm-B73-REFERENCE-NAM-5.0_2Mb.bed", sep = "\t", names = ['chrome','start','end'])

#fold0 = make_site_dict("data/ref/Zm-B73-REFERENCE-NAM-5.0_FOLDDRAFT.txt", 0)
#fold4 = make_site_dict("data/ref/Zm-B73-REFERENCE-NAM-5.0_FOLDDRAFT.txt", 4)

#print("0 FOLD")
#snp_sfs = build_sfsdf(bed_df, sfs_size = 26)
#snp_sfs = file_parse("data/variants/snp_test.vcf", "vcf", bed_df, snp_sfs,  missing_prop = 0/26, heterozygous_prop = 0/26, site_dict = fold0)
#for idx, sfs in snp_sfs.items():
#    bed_row = bed_df.iloc[idx, :]
#    print('\t'.join([str(i) for i in sfs]))

#print("4 FOLD")
#snp_sfs = build_sfsdf(bed_df, sfs_size = 26)
#snp_sfs = file_parse("data/variants/snp_test.vcf", "vcf", bed_df, snp_sfs,  missing_prop = 0/26, heterozygous_prop = 0/26, site_dict = fold4)
#for idx, sfs in snp_sfs.items():
#    bed_row = bed_df.iloc[idx, :]
#    #print(f"{bed_row['chrome']}\t{bed_row['start']}\t{bed_row['end']}\t{' '.join([str(i) for i in sfs])}")
#    print('\t'.join([str(i) for i in sfs]))

#sv_sfs = build_sfsdf(bed_df, sfs_size = 26)
#sv_sfs = file_parse("data/variants/sv_test.txt", "sv", bed_df, sv_sfs,  missing_prop = 0/26, heterozygous_prop = 0/26)
#for idx, sfs in sv_sfs.items():
#    print(idx, sfs)
