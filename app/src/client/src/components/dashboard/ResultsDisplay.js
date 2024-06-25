import React, { useState } from "react";
import {
    Box,
    Button,
    Divider,
    IconButton,
    Paper,
    Typography
} from "@mui/material";
import TypographyTooltip from "../common/TypographyTooltip";
import RefreshIcon from "@mui/icons-material/Refresh";

const ResultsDisplay = ({ 
    results, 
    setResults,
    selectedResultIndex,
    setSelectedResultIndex
}) => {
    const handleSelectResultIndex = (index) => {
        if (selectedResultIndex === index) {
            setSelectedResultIndex(null); // Deselect the result if it is already selected.
        } else {

            setSelectedResultIndex(index); // Select the result.
        }
    };

    const handleRefresh = () => {
        setSelectedResultIndex(null);
        setResults([]);
    };

    return (
        <Paper
            sx={{
                p: 3,
                backgroundColor: "#f0f0f0",
                borderRadius: 2,
                boxShadow: 3,
                width: "100%",
                border: "1px solid #ccc",
                display: "flex",
                flexDirection: "column",
                // minHeight: "385px",
                // maxHeight: "400px",
            }}
        >
            <Box sx={{ flexShrink: 0 }}>
                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="One result resembles a single putative protocluster.">
                            <Typography variant="h6" gutterBottom>
                                Select result
                            </Typography>
                        </TypographyTooltip>
                    </Box>
                    <Box>
                        <IconButton
                            type="button"
                            onClick={handleRefresh}
                            color="primary"
                            sx={{ transform: "translateY(-3px)" }}
                        >
                            <RefreshIcon />
                        </IconButton>
                    </Box>
                </Box>
                <Divider sx={{ mb: 2 }} />
            </Box>
            <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
                <Box sx={{ 
                    flex: "1 1 100%",
                    display: "flex", 
                    // alignItems: "center", 
                    // justifyContent: "center", 
                    padding: "10px", // delete this line when centering component to be implemented
                    backgroundColor: "white", 
                    borderRadius: "4px", 
                    border: "1px solid #bbb",
                    boxShadow: "inset 0 2px 5px rgba(0, 0, 0, 0.1)",
                }}>
                    {selectedResultIndex !== null ? (
                        <Box>
                            {results[selectedResultIndex]["queryType"] === "retrosynthesis" ? (
                                <Box>
                                    <Typography>
                                        Applied rules:
                                    </Typography>
                                    {results[selectedResultIndex]["metaData"]["appliedRules"].map((reaction, index) => (
                                        <Box key={index} sx={{ ml: 2 }}>
                                            <Typography>
                                                {reaction}
                                            </Typography>
                                        </Box>
                                    ))}
                                    <Box sx={{ mt: 1 }}>
                                        <Typography>
                                            {`Number of motifs: ${results[selectedResultIndex]["query"].length}`}
                                        </Typography>
                                    </Box>
                                </Box>
                            ) : (
                                <Box>
                                    <Typography>
                                        {`Category: ${results[selectedResultIndex]["metaData"]["category"]}`}
                                    </Typography>
                                    <Typography>
                                        {`Product: ${results[selectedResultIndex]["metaData"]["product"]}`}
                                    </Typography>
                                    <Typography>
                                        {`Start position: ${results[selectedResultIndex]["metaData"]["startProtocluster"]}`}
                                    </Typography>
                                    <Typography>
                                        {`End position: ${results[selectedResultIndex]["metaData"]["endProtocluster"]}`}
                                    </Typography>
                                    <Typography>
                                        {`Number of motifs: ${results[selectedResultIndex]["query"].length}`}
                                    </Typography>
                                </Box>
                            )}
                        </Box>
                    ) : (
                        <Typography>
                            No result selected.
                        </Typography>
                    )}
                </Box>
                <Box sx={{ 
                    flexGrow: 1, 
                    overflowY: "auto", 
                    // minHeight: "275px",
                    maxHeight: "275px", 
                    flex: "1 1 50%",
                    backgroundColor: "white",
                    boxShadow: "inset 0 2px 5px rgba(0, 0, 0, 0.1)",
                    border: "1px solid #bbb", 
                    borderRadius: "4px",
                    padding: "10px"
                }}>
                    {results.length ? (
                        <Box>
                            {results.map((result, index) => (
                                <Button
                                    key={index}
                                    variant="contained"
                                    color={selectedResultIndex === index ? "secondary" : "primary"}
                                    onClick={() => handleSelectResultIndex(index)}
                                    sx={{ 
                                        width: "100%", 
                                        margin: "5px 0",
                                        display: "flex",
                                        justifyContent: "flex-start",
                                        alignItems: "center",
                                        textAlign: "left",
                                    }}
                                >
                                    <span style={{ marginRight: "10px" }}>{index + 1}</span>
                                    <span>{result["title"]}</span>
                                </Button>
                            ))}
                        </Box>
                    ) : (
                        <Typography>
                            No results available.
                        </Typography>
                    )}
                </Box>
                
            </Box>
        </Paper>
    );
};

export default ResultsDisplay;