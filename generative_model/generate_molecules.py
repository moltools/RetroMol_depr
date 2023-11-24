#!/usr/bin/env python3
import argparse 
import joblib 
import math 
import typing as ty 
from collections import defaultdict 
from random import choice as random_choice
from abc import ABC, abstractmethod

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from retromol.chem import ChemicalReaction, MolecularPattern, mol_to_fingerprint
from retromol.parsing import parse_reaction_rules, parse_molecular_patterns

def cli () -> argparse.Namespace:
    """
    Command line interface.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--model", type=str, required=True, help="Path to the model.")    
    parser.add_argument("-rr", "--reaction-rules", type=str, required=True, help="Path to the reaction rules.")
    parser.add_argument("-cu", "--core-units", type=str, required=True, help="Path to the core units.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output file.")
    return parser.parse_args()

class MCTS:
    def __init__(self, expl_rate: float) -> None:
        self.Q = defaultdict(int) # Total reward 
        self.N = defaultdict(int) # Visits to node
        self.children = {}
        self.expl_rate = expl_rate  # Exploration rate
    
    def choose(self, node: "Node") -> "Node":
        if node.is_terminal(): 
            raise RuntimeError(f"Terminal Node chosen: {node}")
        
        if node not in node.find_children(): 
            return node.find_random_child()
        
        def score(n: Node) -> float:
            if self.N[n] == 0: 
                # Avoid unseen moves.
                return float("-inf")
            return self.Q[n] / self.N[n]
        
        return max(self.children[node], key=score)
    
    def do_rollout(self, node: "Node") -> None:
        path = self._select(node)
        leaf = path[-1]
        self._expand(leaf)
        reward = self._simulate(leaf)
        self._backpropagate(path, reward)
    
    def _select(self, node: "Node") -> ty.List["Node"]:
        path = []
        while True:
            path.append(node)

            if node not in self.children or not self.children[node]: 
                # Node is terminal or unexplored.
                return path
            
            # Subtract known nodes from children node.
            unexplored = self.children[node] - self.children.keys()

            if unexplored:
                n = unexplored.pop()
                path.append(n)
                return path
            
            # Descend one level deeper.
            node = self._ucb_select(node)

    def _expand(self, node: "Node") -> None:
        if node in self.children: 
            return
        
        self.children[node] = node.find_children()

    def _simulate(self, node: "Node") -> None:
        while True:
            if node.is_terminal(): 
                reward = node.reward()
                return reward
            children = node.find_children()
            if not children:
                return node.reward()
            node = node.find_random_child()

    def _backpropagate(self, path: ty.List["Node"], reward: float) -> None:
        for node in reversed(path):
            self.N[node] += 1 
            self.Q[node] += reward 
    
    def _ucb_select(self,node: "Node") -> "Node":
        assert all(n in self.children for n in self.children[node])

        log_N = math.log(self.N[node])

        def ucb(n: Node) -> float:
            return (
                self.Q[n] / self.N[n]
                + self.expl_rate * math.sqrt(log_N / self.N[n])
            )

        return max(self.children[node], key=ucb)

class Node(ABC):
    @abstractmethod
    def find_children(self) -> ty.Set["Node"]:
        ...

    @abstractmethod
    def find_random_child(self) -> ty.Optional["Node"]:
        ...

    @abstractmethod
    def reward(self) -> float:
        ...

    @abstractmethod
    def is_terminal(self) -> bool:
        ...

    @abstractmethod
    def __hash__(self) -> int:
        ...

    @abstractmethod
    def __eq__(self, other: "Node") -> bool:
        ...

class State(Node):
    def __init__(
        self, 
        state: Chem.Mol, 
        prev_state: Chem.Mol,
        depth: int, 
        obj: ty.Callable, 
        reactions: ty.Dict[str, ChemicalReaction],
        motifs: ty.Dict[str, MolecularPattern],
        max_depth: int,
        design_threshold: float,
    ) -> None:
        """
        Constructor for the State class.

        Parameters
        ----------
        state : Chem.Mol
            Current state of the molecule.
        depth : int
            Current depth of the tree.
        obj : ty.Callable
            Objective function.
        reactions : ty.Dict[str, ChemicalReaction]
            Dictionary of reactions.
        motifs : ty.Dict[str, MolecularPattern]
            Dictionary of motifs.
        max_depth : int
            Maximum depth of the tree.
        design_threshold : float
            Threshold for the objective function.
        
        Returns
        -------
        None
        """
        self.state = state
        self.prev_state = prev_state
        self.depth = depth
        self.obj = obj
        self.reactions = reactions
        self.motifs = motifs
        self.max_depth = max_depth
        self.design_threshold = design_threshold
    
    def __repr__(self):
        return f"(State = {self.state}; Depth = {self.depth})"
    
    def find_children(self) -> ty.Set[Node]:
        if self.is_terminal(): 
            return set()
        
        actions = ['macrocyclization', 'add_sugar', 'pks']
        # actions = ['macrocyclization', 'pks']

        child_states = []

        for action in actions:
            if action == "macrocyclization":
                rxn = self.reactions["macrocyclization"].forward_compiled
                results = rxn.RunReactants([self.state])
                action_states = []
                if len(results):
                    for result in results:
                        product = result[0]
                        Chem.SanitizeMol(product)
                        action_states.append(product)
                # Pick best action state.
                if action_states:
                    best = max(action_states, key=lambda x: self.obj([mol_to_fingerprint(x, 4, 2048)]))
                    child_states.append(best)

            elif action == "add_sugar":
                smirks = r"[C:1][OH:2].[OH][C:3]1[C:4]([OH:5])[C:6]([OH:7])[C:8]([OH:9])[C:10]([C:11][OH:12])[O:13]1>>[C:1][O:2][C:3]1[C:4]([OH:5])[C:6]([OH:7])[C:8]([OH:9])[C:10]([C:11][OH:12])[O:13]1"
                sugar_smi = r"C(CO)1OC(O)C(O)C(O)C(O)1"
                sugar = Chem.MolFromSmiles(sugar_smi)
                rxn = rdChemReactions.ReactionFromSmarts(smirks)
                reactants = [self.state, sugar]
                results = rxn.RunReactants(reactants)
                action_states = []
                if len(results):
                    for result in results:
                        product = result[0]
                        Chem.SanitizeMol(product)
                        action_states.append(product)
                # Pick best action state.
                if action_states:
                    best = max(action_states, key=lambda x: self.obj([mol_to_fingerprint(x, 4, 2048)]))
                    child_states.append(best)
            
            else:
                motif_names = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12', 'C1', 'C2', 'C4', 'C6', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D11', 'D12']
                action_states = []
                for motif_name in motif_names:
                    rxn = self.reactions["PKS"].forward_compiled
                    next_motif = self.motifs[motif_name].compiled
                    reactants = [self.state, next_motif]
                    results = rxn.RunReactants(reactants)
                    if len(results):
                        for result in results:
                            product = result[0]
                            Chem.SanitizeMol(product)
                            action_states.append(product)
                # Pick best action state.
                if action_states:
                    best = max(action_states, key=lambda x: self.obj([mol_to_fingerprint(x, 4, 2048)]))
                    child_states.append(best)

        valid_child_states = {
            State(
                state=state,
                prev_state=self.state,
                depth=self.depth + 1,
                obj=self.obj,
                reactions=self.reactions,
                motifs=self.motifs,
                max_depth=self.max_depth,
                design_threshold=self.design_threshold,
            )
            for state in child_states
        }

        return valid_child_states

    def find_random_child(self) -> ty.Optional["State"]:
        if self.is_terminal():
            return None
        return random_choice(list(self.find_children()))

    def reward(self) -> float:
        # fp_prev = mol_to_fingerprint(self.prev_state, 4, 2048)
        fp_curr = mol_to_fingerprint(self.state, 4, 2048)
        # return_val = self.obj([fp_curr]) - self.obj([fp_prev])
        # print(return_val)
        # return return_val
        return self.obj([fp_curr])

    def is_terminal(self) -> bool:
        fp = mol_to_fingerprint(self.state, 4, 2048)

        if self.depth >= self.max_depth:
            return True
        elif self.obj([fp]) >= self.design_threshold:
            return True
        # elif len(self.find_children()) == 0:
        #     return True
        else:
            return False

    def __eq__(self, other: "State") -> bool:
        return self.__hash__() == other.__hash__()
    
    def __hash__(self):
        fp = mol_to_fingerprint(self.state, 4, 2048)
        fp = np.append(fp, self.depth)
        return hash(tuple(fp))

def main () -> None:
    """
    Main function.
    """
    args = cli()
    model = joblib.load(args.model)
    x0 = Chem.MolFromSmiles(r"[C][C][C](=[O])[OH]")
    Chem.SanitizeMol(x0)

    def func(x: np.array) -> float:
        return model.predict_proba(x)[0][1]
    
    reactions = parse_reaction_rules(args.reaction_rules)
    reactions_dict = {reaction.name: reaction for reaction in reactions if reaction.name in ["PKS", "macrocyclization"]}
    motifs = parse_molecular_patterns(args.core_units, as_smiles=True)
    def is_pks(name: str) -> bool:
        return name[0] in ["A", "B", "C", "D"] and name[1:].isdigit()
    motifs_dict = {motif.name: motif for motif in motifs if is_pks(motif.name)}
    
    process = MCTS(expl_rate=1.0)
    state = State(x0, None, 0, func, reactions_dict, motifs_dict, max_depth=10, design_threshold=0.5)

    solution = [state.state]
    while True:
        for i in range(10):
            process.do_rollout(state)
            print(f"{i}".zfill(3).format(i), end="\r")

        state = process.choose(state)

        state_fp = mol_to_fingerprint(state.state, 4, 2048)
        state_pred = func([state_fp])
        print(state)
        print(Chem.MolToSmiles(state.state), state_pred)

        solution.append(state.state)
        if state.is_terminal(): 
            break   
        children = state.find_children()
        if not children:
            break

    print(len(solution))
    if solution[-1] is not None:
        print(Chem.MolToSmiles(solution[-1]))
        fp = mol_to_fingerprint(solution[-1], 4, 2048)
        print(func([fp]))
    else:
        print("No solution found.")

    exit(0)

if __name__ == "__main__":
    main()