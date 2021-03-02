import os
import glob
import pandas as pd
import argparse


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", required=True,
                    help="Enter into the folder")
    parser.add_argument("-omega", "--omega_path", required=True,
                    help="Enter the OMEGA path")
    parser.add_argument("-rocs", "--rocs_path", required=True,
                    help="Enter the ROCS path")
    parser.add_argument("-temp_lig", "--template_ligand_path", required=True,
                    help="Enter the template ligand PDB path")
    parser.add_argument("-mol2params", required=True,
                    help="Enter the template ligand PDB path")    
    return parser.parse_args()



def conf_gen(smi, maxconfs=2000):

    smi_prefix = os.path.splitext(os.path.basename(smi))[0]

    cmd = f'{OMEGA} -in {smi} -out OMEGA/{smi_prefix}_omega.sdf \
                -prefix OMEGA/{smi_prefix}_omega -warts true \
                -maxconfs {maxconfs} -strict false'

    os.system(cmd)


# ligand alignment using ROCS openeye

def lig_alignment(conformer, template_database, rocs_maxconfs_output=100):

    sdf_prefix = os.path.basename(os.path.splitext(conformer)[0]).split('_')[0]
    
    for template in template_database:
        template_id = "_".join(os.path.basename(template).split("_")[0:3])

        cmd = f'{ROCS} -dbase {conformer} -query {template} \
               -prefix {sdf_prefix}_{template_id}_rocs -oformat mol2 \
               -maxconfs {rocs_maxconfs_output} -outputquery false \
               -qconflabel title -outputdir ROCS/'

        os.system(cmd)



# Combine each report file into single file


def combine_report_files(report_file):

    data = []

    for rpt in report_file:
        target_template_name = os.path.basename(rpt).replace("_1.rpt", "")

        rpt = pd.read_csv(rpt, sep='\t')
        rpt = rpt.loc[:, ~rpt.columns.str.match('Unnamed')]

        rpt["ShapeQuery"] = target_template_name
        
        data.append(rpt)
    
    data = pd.concat(data)
    data = data.sort_values(by=['TanimotoCombo'], ascending=False)
    data["Rank"] = range(1, data.shape[0]+1)

    data.to_csv('ROCS/single_report_file_sorted.csv', index=False)
    
    data_100 = data.iloc[:100, :][["Name", "ShapeQuery"]]
    data_100.to_csv("ROCS/top_100.txt", index=False)
    

# Seperate the top 100 conformer hits from ROCS alignment sdf files


def sep_hits_from_rocs_sdf_file(top_100_hits_txt_path):

    def extract_mol2_conf(file, conf_labels):
        with open(f"ROCS/{file}_hits_1.mol2", "r") as f:
            f = f.read().split("@<TRIPOS>MOLECULE\n")
            f = {v.split("\n")[0]:"@<TRIPOS>MOLECULE\n" + v for v in f}
            
            for label in conf_labels:
                with open(f"top_100_conf/{label}_{file}_hits.mol2", "w") as fwr:
                    fwr.write(f[label])
                    fwr.close()

    data = pd.read_csv(top_100_hits_txt_path)
    data = data.groupby("ShapeQuery")["Name"].apply(list).to_dict()

    for template, conformers in data.items():
        extract_mol2_conf(template, conformers)


# Convert SDF file to PDB/PARAMS for Rosetta input


def sdftomol2(mol2params, top_hits_sdf_path):

    for file in top_hits_sdf_path:

        prefix = file.split("/")[-1].split(".")[0]
        cmd = f'python {mol2params} -s {file} --prefix=mol2params/{prefix}'

        os.system(cmd)



if __name__ == "__main__":

    args = args()

    OMEGA = args.omega_path
    ROCS = args.rocs_path
    MOL2PARAMS = args.mol2params

    template_lig_library = args.template_ligand_path
    template_lig_library = glob.glob(f"{template_lig_library}/*.pdb")

    os.chdir(args.folder)

    os.mkdir("OMEGA")
    os.mkdir("ROCS")
    os.mkdir("top_100_conf")
    os.mkdir("mol2params")

    smiles = glob.glob("*.smi")[0]
    conf_gen(smiles, maxconfs=1000)

    sdf = glob.glob("OMEGA/*.sdf")[0]
    lig_alignment(sdf, template_lig_library, rocs_maxconfs_output=30)

    combine_report_files(glob.glob("ROCS/*.rpt"))
    sep_hits_from_rocs_sdf_file("ROCS/top_100.txt")

    sdftomol2(MOL2PARAMS, glob.glob("top_100_conf/*.mol2"))
