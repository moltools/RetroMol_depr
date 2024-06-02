import React from "react";
import { Link } from "react-router-dom";

import Page from "../components/common/Page";

const CineMol = () => {
    return (
        <Page>
            <h1>Under construction</h1>
            <p>This page is under construction. Please check back later.</p>
            <p>Click <Link to="/">here</Link> to return to the homepage.</p>
        </Page>
    );
};

export default CineMol;