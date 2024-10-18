import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import { Parser, SvgDrawer } from 'smiles-drawer';
import { Canvas } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';

// const smiles = "CC(=O)Oc1ccccc1C(=O)O";
// const smiles = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C";
const smiles = "CCCCCCCCCC(=O)NC(CC1=CNC2=CC=CC=C21)C(=O)NC(CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC3C(OC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C";

const atomColors = {
    "H": "#FFFFFF",
    "C": "#DFDFDF",
    "N": "#0000FF",
    "O": "#FF0000",
    "F": "#00FF00",
    "Cl": "#00FF00",
    "Br": "#A52A2A",
    "I": "#9400D3",
    "P": "#FFA500",
    "S": "#FFD700",
};

function getAtomColor(atomSymbol) {
    const defaultColor = "#808080";

    if (atomSymbol in atomColors) {
        return atomColors[atomSymbol];
    } else if (atomSymbol.toUpperCase() in atomColors) {
        return atomColors[atomSymbol.toUpperCase()];
    };

    return defaultColor;
};

function getMidpoint(source, target) {
    return [(source[0] + target[0]) / 2, (source[1] + target[1]) / 2, 0];
};

const Atom = ({ position, color }) => {
    return (
        <mesh 
            position={[position[0], position[1], 0]}
        >
            <circleGeometry args={[6, 32]} />
            <meshStandardMaterial color={color} />
        </mesh>
    );
};
  
const Bond = ({ sourcePosition, targetPosition, color }) => {
    const getAngle = (source, target) => {
        return Math.atan2(target[1] - source[1], target[0] - source[0]) + Math.PI / 2;
    };

    const getLength = (source, target) => {
        return Math.sqrt(Math.pow(target[0] - source[0], 2) + Math.pow(target[1] - source[1], 2));
    };

    return (
        <mesh 
            position={getMidpoint(sourcePosition, targetPosition)}
            rotation={[0, 0, getAngle(sourcePosition, targetPosition)]}
        >
            <cylinderGeometry args={[6, 6, getLength(sourcePosition, targetPosition), 16]} />
            <meshStandardMaterial color={color} />
        </mesh>
    );
};

const Molecule = ({ atoms, bonds }) => {
    const [zoom, setZoom] = useState(100);

    useEffect(() => {
        // Calculate bounds.
        let minX = Infinity;
        let minY = Infinity;
        let maxX = -Infinity;
        let maxY = -Infinity;

        atoms.forEach((atom) => {
            minX = Math.min(minX, atom.position[0]);
            minY = Math.min(minY, atom.position[1]);
            maxX = Math.max(maxX, atom.position[0]);
            maxY = Math.max(maxY, atom.position[1]);
        }, []);

        // Calculate zoom.
        const width = (maxX - minX) * 1.25;
        const height = (maxY - minY) * 1.25;
        const zoom = (1 / (Math.max(width, height) / 2) * 100);

        setZoom(zoom);
    }, [atoms]);

    return (
        <Canvas 
            orthographic
            camera={{ 
                position: [0, 0, 100], 
                zoom: zoom 
            }}
            style={{ 
                width: "100%", 
                height: "100%",
                backgroundColor: "#f0f0f0" 
            }}
        >
            <ambientLight />
            <OrbitControls enableRotate={false} />
            {bonds.map((bond, index) => {
                let source = atoms[bond.source];
                let target = atoms[bond.target];
                let midpointPosition = getMidpoint(source.position, target.position);
                
                return (
                    <Bond 
                        key={index} 
                        sourcePosition={source.position}
                        targetPosition={midpointPosition}
                        color={source.color}
                    />
                );
            })}
            {atoms.map((atom, index) => (
                <Atom 
                    key={index} 
                    position={atom.position} 
                    color={atom.color}
                />
            ))}
        </Canvas>
    );
};

const Home = () => {
    const [atoms, setAtoms] = useState([]);
    const [bonds, setBonds] = useState([]);

    useEffect(() => {
        let parseTree = Parser.parse(smiles);
        const drawer = new SvgDrawer();
        const preprocessor = drawer.preprocessor;
        preprocessor.initDraw(parseTree);
        preprocessor.processGraph();

        let atoms = [];
        let bonds = [];
        let centroid = { x: 0, y: 0 };

        preprocessor.graph.vertices.forEach((vertex) => {
            const xPos = vertex.position.x;
            const yPos = vertex.position.y;
            const atomSymbol = vertex.value.element;
            atoms.push({ 
                "position": [xPos, yPos, 0.0],
                "color": getAtomColor(atomSymbol) 
            });
            centroid.x += xPos;
            centroid.y += yPos;
        });

        atoms.forEach((atom) => {
            atom.position[0] -= centroid.x / atoms.length;
            atom.position[1] -= centroid.y / atoms.length;
        });
        setAtoms(atoms);
        
        preprocessor.graph.edges.forEach((edge) => {
            const sourceId = edge.sourceId;
            const targetId = edge.targetId;
            bonds.push({ "source": sourceId, "target": targetId });
            bonds.push({ "source": targetId, "target": sourceId });
        });
        setBonds(bonds);

    }, []);

    return (
        <div style={{ display: "flex", flexDirection: "column" }}>
            {/* <h1>RetroMol</h1>
            <Link to="/">back</Link>
            <Link to="/retromol/workspace">workspace</Link>
            <div style={{ width: "200px", height: "200px" }}>
                <Molecule 
                    atoms={atoms} 
                    bonds={bonds} 
                />
            </div> */}
            RetroMol is currently being prepared for submission. Please check back later.
        </div>
    );
};

export default Home;