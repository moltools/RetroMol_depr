import React from "react";
import { Container, CssBaseline, Box } from "@mui/material";

const Page = ({ children, ...props }) => {
    return (
        <Container {...props}>
            <CssBaseline />
            <Box my={4}>
                {children}
            </Box>
        </Container>
    );
};

export default Page;