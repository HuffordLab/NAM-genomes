#Manually declage the GFF Star indexes as well as the prefixes that we want to use
declare -a FASTA=("/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/B73_NAMassembly/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/Ki11_NAMassembly/Zm-Ki11-REFERENCE-NAM-1.0/Zm-Ki11-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/P39_NAMassembly/Zm-P39-REFERENCE-NAM-1.0/Zm-P39-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/NC350_NAMassembly/Zm-NC350-REFERENCE-NAM-1.0/Zm-NC350-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/B97_NAMassembly/Zm-B97-REFERENCE-NAM-1.0/Zm-B97-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/CML103_NAMassembly/Zm-CML103-REFERENCE-NAM-1.0/Zm-CML103-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/CML228_NAMassembly/Zm-CML228-REFERENCE-NAM-1.0/Zm-CML228-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/CML247_NAMassembly/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/CML277_NAMassembly/Zm-CML277-REFERENCE-NAM-1.0/Zm-CML277-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/CML322_NAMassembly/Zm-CML322-REFERENCE-NAM-1.0/Zm-CML322-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/CML333_NAMassembly/Zm-CML333-REFERENCE-NAM-1.0/Zm-CML333-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/CML52_NAMassembly/Zm-CML52-REFERENCE-NAM-1.0/Zm-CML52-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/CML69_NAMassembly/Zm-CML69-REFERENCE-NAM-1.0/Zm-CML69-REFERENCE-NAM-1.0.fasta"           
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/HP301_NAMassembly/Zm-HP301-REFERENCE-NAM-1.0/Zm-HP301-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/Il14H_NAMassembly/Zm-Il14H-REFERENCE-NAM-1.0/Zm-Il14H-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/Ki3_NAMassembly/Zm-Ki3-REFERENCE-NAM-1.0/Zm-Ki3-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/M162W_NAMassembly/Zm-M162W-REFERENCE-NAM-1.0/Zm-M162W-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/M37W_NAMassembly/Zm-M37W-REFERENCE-NAM-1.0/Zm-M37W-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/Mo18W_NAMassembly/Zm-Mo18W-REFERENCE-NAM-1.0/Zm-Mo18W-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/NC358_NAMassembly/Zm-NC358-REFERENCE-NAM-1.0/Zm-NC358-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/Oh7B_NAMassembly/Zm-Oh7B-REFERENCE-NAM-1.0/Zm-Oh7B-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/Oh43_NAMassembly/Zm-Oh43-REFERENCE-NAM-1.0/Zm-Oh43-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/Ky21_NAMassembly/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/Tx303_NAMassembly/Zm-Tx303-REFERENCE-NAM-1.0/Zm-Tx303-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/Tzi8_NAMassembly/Zm-Tzi8-REFERENCE-NAM-1.0/Zm-Tzi8-REFERENCE-NAM-1.0.fasta"
                  "/panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/Ms71_NAMassembly/Zm-Ms71-REFERENCE-NAM-1.0/Zm-Ms71-REFERENCE-NAM-1.0.fasta")




declare -a QUERY=("B73"
                  "Ki11"
                  "P39"
                  "NC350"
                  "B97"
                  "CML103"
                  "CML228"
                  "CML247"
                  "CML277"
                  "CML322"
                  "CML333"
                  "CML52"
                  "CML69"
                  "HP301"
                  "Il14H"
                  "Ki3"
                  "M162W"
                  "M37W"
                  "Mo18W"
                  "NC358"
                  "Oh7B"
                  "Oh43"
                  "Ky21"
                  "Tx303"
                  "Tzi8"
                  "Ms71")

for j in {0..25}; do
    for i in {0..25}; do
      if [ "$i" -ne "$j" ]; then
        (echo nucmer --mum -c 1000 -p ${QUERY[$j]}_${QUERY[$i]}_c250 ${FASTA[$j]} ${FASTA[$i]})
      fi
    done
done
