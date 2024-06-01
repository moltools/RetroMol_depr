import React from "react";
import {
    Box,
    Divider,
    Paper,
    Typography
} from "@mui/material";
import TypographyTooltip from "../common/TypographyTooltip";

const ResultsDisplay = ({ results }) => {
    return (
        <Paper
            sx={{
                p: 3,
                backgroundColor: "#f0f0f0",
                borderRadius: 2,
                boxShadow: 3,
                width: "100%",
                border: "1px solid #ccc"
            }}
        >
            <Box sx={{ display: "flex", gap: 1 }}>
                <TypographyTooltip title="One result resembles a single putative protocluster.">
                    <Typography variant="h6" gutterBottom>
                        Select result
                    </Typography>
                </TypographyTooltip>
            </Box>
            <Divider sx={{ mb: 2 }} />
            {results.length ? (
                <Box>
                    {results.map((result, index) => (
                        <Box key={index}>
                            <Typography>
                                {JSON.stringify(result, null, 2)}
                            </Typography>
                        </Box>
                    ))}
                </Box>
            ) : (
                <Typography>
                No results available.
                </Typography>
            )}
        </Paper>
    );
};

export default ResultsDisplay;