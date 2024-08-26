**Reward Function for Drug Discovery Reinforcement Learning via Symbolic Feedback (RLSF)**

Input: Array of reactants, product. They must be written in SMILES format.
Output: Binary reward vector representing precise location of errors in molecule

Sample Usage: Valid Reaction
Input:
- product_smiles = 'CC(=O)OCC1=CC=CC=C1C(=O)O'
- reactant_smiles_array = ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CCO']
Output:
- The predicted reactants are valid and produce the correct product.

Sample Usage: Invalid Reaction
Input:
- product_smiles = 'CC(=O)OCC1=CC=CC=C1C(=O)O'
- reactant_smiles_array = ['CC(=O)O[C1]=CC=CC=C1C(=O)O', 'CC1=CC=CC1=O']
Output:
- The reactants do not combine to produce the given product.
- Product Error Vector: [1, 1, 1, 1, 1, 1, 1]
- Reactant 1 Error Vector: [0, 1, 1, 1, 1, 1, 0]
- Reactant 2 Error Vector: [0, 1, 1, 1, 1, 0, 1]
