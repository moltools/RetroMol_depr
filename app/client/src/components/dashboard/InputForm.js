import React, { useState } from "react";
import { 
    Button, 
    Box, 
    Divider,
    FormControlLabel, 
    IconButton,
    Paper, 
    Radio, 
    RadioGroup, 
    TextField, 
    Typography 
} from '@mui/material'
import TypographyTooltip from "../common/TypographyTooltip";
import RefreshIcon from '@mui/icons-material/Refresh';

const smilesErythromycin = "CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O";

const AntismashLogo = ({ style }) => {
    // Make the logo stand out from the background.
    const defaultStyle = { boxShadow: "0 0 10px 5px white" };
  
    return (
        <a href="https://antismash.secondarymetabolites.org/" target="_blank" rel="noreferrer">
            <img 
                src={`${process.env.PUBLIC_URL + "/imgs/antismash_logo.svg"}`} 
                alt="Antismash logo" 
                style={{ ...defaultStyle, ...style }}
            />
        </a>
    );
};

const InputForm = ({ onSubmit }) => {
    const jsonFileNameDefault = "Upload JSON file";

    const [selectedInputType, setSelectedInputType] = useState("smiles");
    const [inputValue, setInputValue] = useState("");
    const [jsonContent, setJsonContent] = useState("");
    const [jsonFileName, setJsonFileName] = useState(jsonFileNameDefault);

    const handleRefresh = () => {
        setInputValue("");
        setJsonContent("");
        setJsonFileName(jsonFileNameDefault);
        setSelectedInputType("smiles");
    };

    const handleInputTypeChange = (event) => {
        setInputValue("");
        setJsonContent("");
        setJsonFileName(jsonFileNameDefault);
        setSelectedInputType(event.target.value);
    };

    const handleInputChange = (event) => {
        setInputValue(event.target.value);
    };

    const handleFileChange = (event) => {
        const file = event.target.files[0];
        setJsonFileName(file.name);
        const reader = new FileReader();
        reader.onload = (event) => {
            setJsonContent(event.target.result);
        };
        reader.readAsText(file);
    };

    const handleFormSubmit = (event) => {
        event.preventDefault();
        if (selectedInputType === "json") {
            onSubmit({
                inputType: selectedInputType,
                inputValue: jsonContent
            });
        } else {
            onSubmit({
                inputType: selectedInputType,
                inputValue: inputValue
            });
        };
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
                minHeight: "385px",
                maxHeight: "385px",
            }}
        >
            <form onSubmit={handleFormSubmit}>
                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="AntiSMASH results can be parsed by either passing along a job ID or uploading a result in JSON format.">
                            <Typography variant="h6" gutterBottom>
                                Select input type
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
                <RadioGroup
                    value={selectedInputType}
                    onChange={(event) => handleInputTypeChange(event)}
                    sx={{ mb: 2 }}
                >
                    <FormControlLabel 
                        value="smiles" 
                        control={<Radio />} 
                        label="SMILES" 
                        sx={{ marginLeft: "25px" }}
                    />
                    <Box>
                        <AntismashLogo style={{ width: "25px", transform: "translateY(7px)" }} />
                        <FormControlLabel 
                            value="jobId" 
                            control={<Radio />} 
                            label="Job ID" 
                            sx={{ marginLeft: "0px" }}
                        />
                    </Box>
                    <Box>
                        <AntismashLogo style={{ width: "25px", transform: "translateY(7px)" }} />
                        <FormControlLabel 
                            value="json" 
                            control={<Radio />} 
                            label="JSON" 
                            sx={{ marginLeft: "0px" }}
                        />
                    </Box>
                </RadioGroup>
                {selectedInputType !== "json" && (
                    <TextField
                        sx={{ width: "100%" }}
                        label={selectedInputType === "smiles" ? "Enter SMILES string" : "Enter Job ID"}
                        value={inputValue}
                        onChange={handleInputChange}
                        margin="normal"
                        variant="outlined"
                    />
                )}
                {selectedInputType === "json" &&(
                    <Button
                        variant="outlined"
                        component="label"
                        sx={{ 
                            mt: 3, 
                            width: "100%", 
                            height: "56px", // Align JSON file upload button with text field.
                            transform: "translateY(-8px)" // Align JSON file upload button with text field.
                        }}
                        disabled={jsonContent !== ""}
                    >
                        {jsonFileName}
                        <input
                            type="file"
                            hidden
                            accept=".json"
                            onChange={handleFileChange}
                        />
                    </Button>
                )}
                <Box sx={{ display: "flex", justifyContent: "flex-end", gap: 2 }}>
                    <Button
                        type="button"
                        variant="contained"
                        sx={{ mt: 2, width: "150px" }}
                        color="primary"
                        onClick={() => {
                            setSelectedInputType("smiles");
                            setInputValue(smilesErythromycin);
                        }}
                    >
                        Example
                    </Button>
                    <Button 
                        type="submit" 
                        variant="contained" 
                        sx={{ mt: 2, width: "150px" }}>
                        Submit
                    </Button>
                </Box>
            </form>
        </Paper>
    );
};

export default InputForm;