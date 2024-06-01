import React from "react";
import { Link } from "react-router-dom";

import Page from "../components/common/Page";

const NotFoundPage = () => {
    return (
        <Page>
            <h1>Page not found</h1>
            <p>The page you are looking for does not exist.</p>
            <p>Click <Link to="/">here</Link> to return to the homepage.</p>
        </Page>
    );
};

export default NotFoundPage;