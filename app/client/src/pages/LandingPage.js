import React from "react";
import { Link } from "react-router-dom";

import Page from "../components/common/Page";

const LandingPage = () => {
    return (
        <Page>
            <h1>MolTools</h1>
            <p>Welcome to MolTools, a collection of tools for cheminformatics and molecular modeling of natural.</p>
            <p><Link to="/cinemol">CineMol</Link>: a direct-to-SVG 3D small molecule drawer.</p>
            <p><Link to="/retromol">RetroMol</Link>: a web-based retrosynthesis tool for modular natural products.</p>
        </Page>
    );
};

export default LandingPage;