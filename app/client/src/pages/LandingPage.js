import React from "react";
import { Link } from "react-router-dom";

import Page from "../components/common/Page";

const LandingPage = () => {
    return (
        <Page>
            <h1>RetroMol</h1>
            <p>Welcome to RetroMol, a web-based retrosynthesis tool for modular natural products.</p>
            <p>Click <Link to="/dashboard">here</Link> to access the dashboard.</p>
        </Page>
    );
};

export default LandingPage;