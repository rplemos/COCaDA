<p align="center">
  <img src="https://github.com/user-attachments/assets/57f17d48-baf4-4bed-97ab-6b817e31dc26" alt="COCaDA_logomenor">
</p>

# COCαDA - Large Scale Protein Interatomic Contact Optimization by Cα Distance Matrices

---

## 🔬 Description

COCαDA (Contact Optimization by alpha-Carbon Distance Analysis) optimizes the calculation of atomic interactions in proteins, by using a set of fine-tuned **Cα distances** between every pair of aminoacid residues.
The code includes a customized parser for both **PDB** and **CIF** files, with support for large files, residue and atom filtering, and geometric analysis (e.g., centroids and normal vectors for aromatic residues). Users can also define their own contact distance cutoffs via the [`contact_distances.json`](contact_distances.json) configuration file.

### 🔍 Contact types detected:
- Hydrophobic
- Hydrogen Bond
- Attractive
- Repulsive
- Disulfide Bond
- Salt Bridge
- Aromatic Stacking

---

## 🚀 Features

- ⚡ **Fast Processing**: COCaDA averages 2.5x faster processing times against Fixed Cutoffs definitions, and 6x faster against Biopython's `NeighborSearch`.
- 📂 **PDB and CIF Parsing**: Efficient parsing for both PDB and CIF files to extract atomic and residue information.
- 🧼 **Residue Filtering**: Ignores low-quality atoms, water molecules, and problematic/erroneous information.
- 🔬 **Interaction Analysis**: Identifies contacts between specified types of atoms based on predefined conditions.
- ⚙️ **User-defined Distance Cutoffs**: Predefined distance cutoffs for all seven contact types can be easily changed by the user.
- 🌀 **Aromatic Stacking Detection**: Computes centroids and normal vectors for aromatic residues, to determine aromatic stacking contacts.
- 🧠 **Multi-Core Processing**: Parallel batch processing across any combination of CPU cores.
- 📊 **CSV Output**: Clean, structured results ideal for post-analysis and exploratory tasks.

---

## 📦 Installation

### 🔧 Prerequisites

- Python ≥ 3.x
- `numpy` ≥ 2.0.1
- `psutil` ≥ 6.0.0 *(required only for multi-core mode)*

### 📥 Setup

1. Clone the repository:
   ```sh
   git clone https://github.com/LBS-UFMG/COCaDA.git
   ```

2. Navigate into the project directory:
   ```sh
   cd COCaDA
   ```

3. Set up a virtual environment (optional but recommended):
    ```sh
    python -m venv venv
    source venv/bin/activate
    ```

4. Install dependencies:
    ```sh
    pip install -r requirements.txt
    ```
    
---

## 🛠️ Usage
### ▶️ To run COCaDA:

1. Ensure you are in the project directory and the virtual environment is activated (if used).

2. Run COCaDA with the path to your folder or file:
    ```sh
    python cocada.py -f path_to_files/*.cif [-m [CORES]] [-o [OUTPUT_DIR]] [-d] [-h]
    ```
    
### ⚙️ Parameters

| Flag | Description |
|------|-------------|
| `-f`, `--files` | **(Required)** Path(s) to `.cif` or `.pdb` file(s). Wildcards (e.g., `*.cif`) are accepted. |
| `-m`, `--mode` | **(Optional)** Enables Multi-Core mode. <br>• No value = uses all available cores. <br>• `-m 2` = use core 2 only. <br>• `-m 1-4` = use cores 1 to 4. <br>• `-m 0,2,5` = use specified cores. |
| `-o`, `--output` | **(Optional)** Outputs detailed results to CSV-formatted files. <br>• No value = saves to `./outputs`. <br>• `-o custom_folder` = saves to specified folder. |
| `-d`, `--distances` | **(Optional)** Reads custom cutoff distances from the [`contact_distances.json`](contact_distances.json) configuration file. <br>By default, predefined values are used. |
| `-h`, `--help` | Shows help and usage instructions. |

<p align="center">
  <img src="https://github.com/user-attachments/assets/e6ee3153-9b2c-47a2-b57c-4f63f206a34b" alt="COCaDA_overview">
</p>

---

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🧾 Citation

**LEMOS, Rafael P.; MARIANO, Diego; SILVEIRA, Sabrina A.; MELO-MINARDI, Raquel C. de.**  
*COCαDA - Large-Scale Protein Interatomic Contact Cutoff Optimization by Cα Distance Matrices.*  
Proceedings of the XVII Brazilian Symposium on Bioinformatics (BSB), 17, pp. 59–70, 2024.  
DOI: [https://doi.org/10.5753/bsb.2024.245545](https://doi.org/10.5753/bsb.2024.245545)

---

## 👤 Contact
For any questions or issues, please contact:

**Rafael Pereira Lemos**  
PhD Student in Bioinformatics  
Federal University of Minas Gerais  

- 📧 Email: [rafaellemos42@gmail.com](mailto:rafaellemos42@gmail.com)  
- 🔗 GitHub: [@rplemos](https://github.com/rplemos)

---

## 🧠 Contributions and Acknowledgements
 - Prof. Raquel Cardoso de Melo Minardi, UFMG;
 - Prof. Sabrina de Azevedo Silveira, UFV;
 - Dr. Diego César Batista Mariano, UFMG;
 - All the 'Laboratory of Bioinformatics and Systems' team.
