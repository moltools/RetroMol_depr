import React from "react";
import { Box } from "@mui/material";

const Status = ({ statusName, status }) => (
    <Box
        variant="contained"
        sx={{
            mb: 1,
            backgroundColor: status ? "#78c419" : "#ec462e",
            transition: "background-color 0.5s ease",
            color: "white",
            textAlign: "center",
            padding: "8px",
            borderRadius: "4px",
            cursor: "default",
            userSelect: "none",
            color: "#222222",
            mx: 1
        }}
    >
        {status ? `${statusName}: Online` : `${statusName}: Offline`}
    </Box>
);

export default Status;