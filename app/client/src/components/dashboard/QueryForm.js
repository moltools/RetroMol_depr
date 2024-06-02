import React, { useEffect, useState } from "react";
import { Box, Button, Divider, FormControl, FormControlLabel, IconButton, InputLabel, MenuItem, MenuProps, Paper, 
    PaperProps, Select, Switch, TextField, Tooltip, Typography ,
    Radiogroup, Radio, RadioGroup
} from "@mui/material";
import { DragDropContext, Droppable, Draggable } from "react-beautiful-dnd";
import TypographyTooltip from "../common/TypographyTooltip";
import RefreshIcon from "@mui/icons-material/Refresh";
import CloseIcon from "@mui/icons-material/Close";
import AddIcon from "@mui/icons-material/Add";
import SearchIcon from "@mui/icons-material/Search";
import { toast } from "react-toastify";
import MultiSelect from "../common/MultiSelect";

const defaultMotif = {
    motifType: "polyketide",
    polyketideType: "Any",
    polyketideDecor: "Any",
    peptideSource: "Any",
    peptideCid: "Any",
};

const QueryForm = ({ results, selectedResultIndex, columns, setColumns, submit }) => {
    
    // query settings
    const [queryType, setQueryType] = useState("match"); // match or query
    const [alignmentType, setAlignmentType] = useState("global"); // global or local
    const [gapPenalty, setGapPenalty] = useState(2);
    const [endGapPenalty, setEndGapPenalty] = useState(1);
    const [queryHasLeadingModules, setQueryHasLeadingModules] = useState(false);
    const [queryHasTrailingModules, setQueryHasTrailingModules] = useState(false);
    
    const [allBioactivityLabels, setAllBioactivityLabels] = useState([]);
    const [selectedBioactivityLabels, setSelectedBioactivityLabels] = useState([]);
    const [allOrganismLabels, setAllOrganismLabels] = useState([]);
    const [selectedOrganismLabels, setSelectedOrganismLabels] = useState([]);

    const [maxNumMatches, setMaxNumMatches] = useState(50);
    const [queryAgainstMolecules, setQueryAgainstMolecules] = useState(true);
    const [queryAgainstProtoclusters, setQueryAgainstProtoclusters] = useState(false);
    const [minMatchLength, setMinMatchLength] = useState(1);
    const [maxMatchLength, setMaxMatchLength] = useState(100);

    const handleRefresh = () => {
        setQueryType("match");
        setAlignmentType("global");
        setGapPenalty(2);
        setEndGapPenalty(1);
        setQueryHasLeadingModules(false);
        setQueryHasTrailingModules(false);
        setSelectedBioactivityLabels([]);
        setSelectedOrganismLabels([]);
        setMaxNumMatches(50);
        setQueryAgainstMolecules(true);
        setQueryAgainstProtoclusters(false);
        setMinMatchLength(1);
        setMaxMatchLength(100);
        
        // when you start querying right away, the selectedResultIndex is null
        if (selectedResultIndex === null) {
            setColumns([]);
            return;
        };
        setColumns(results[selectedResultIndex]);
    };

    const onDragEnd = (result) => {
        const { destination, source } = result;
        if (!destination) {
            return;
        };
        const newColumns = Array.from(columns);
        const [movedItem] = newColumns.splice(source.index, 1);
        newColumns.splice(destination.index, 0, movedItem);
        setColumns(newColumns);
    };

    const handleDeleteColumn = (index) => {
        const newColumns = columns.filter((_, i) => i !== index);
        setColumns(newColumns);
    };


    const handleAddColumn = () => {
        const newColumns = [...columns, [defaultMotif]];
        setColumns(newColumns);
    };

    const addMotifToColumn = (index) => () => {
        const newColumns = Array.from(columns);
        newColumns[index] = [...newColumns[index], defaultMotif];
        setColumns(newColumns);
    };

    const fetchBioactivityLabels = async () => {
        try {
            const response = await fetch("/api/fetch_bioactivity_labels", {
                method: "GET",
                headers: { "Content-Type": "application/json" }
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();

            if (json.status === "success") {
                setAllBioactivityLabels(json.payload.bioactivityLabels);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
        } catch (error) {
            toast.error(error.message);
        };
    };

    const fetchOrganismLabels = async () => {
        try {
            const response = await fetch("/api/fetch_organism_labels", {
                method: "GET",
                headers: { "Content-Type": "application/json" }
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();

            if (json.status === "success") {
                setAllOrganismLabels(json.payload.organismLabels);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
        } catch (error) {
            toast.error(error.message);
        };
    };

    useEffect(() => {
        fetchBioactivityLabels();
        fetchOrganismLabels();
        if (selectedResultIndex !== null) {
            const selectedResult = results[selectedResultIndex];
            setColumns(selectedResult);
        } else {
            setColumns([]);
        };
    }, [selectedResultIndex]);

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
            }}
        >
            <Box sx={{ flexShrink: 0 }}>
                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="The query is used to search for similar protoclusters through pairwise matching or pattern matching.">
                            <Typography variant="h6" gutterBottom>
                                Build query
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
                <Box>
                    <Box 
                        style={{
                            minHeight: "50px", 
                            backgroundColor: "#555",
                            borderTopLeftRadius: "4px",
                            borderTopRightRadius: "4px",
                        }}
                    >
                        <Button
                            variant="contained"
                            color="primary"
                            sx={{ margin: "8px", width: "150px" }}
                            onClick={handleAddColumn}
                        >
                            Add motif
                        </Button>
                    </Box>
                    <Box sx={{ 
                        flexGrow: 1, 
                        overflowY: "auto", 
                        flex: "1 1 50%",
                        minHeight: "275px",
                        backgroundColor: "white",
                        boxShadow: "inset 0 2px 5px rgba(0, 0, 0, 0.1)",
                        border: "1px solid #bbb", 
                        borderTopLeftRadius: "0",
                        borderTopRightRadius: "0",
                        borderBottomLeftRadius: "4px",
                        borderBottomRightRadius: "4px",
                        pt: 2,
                    }}>
                        {(columns && columns.length) ? (
                            <Box 
                                sx={{ 
                                    overflowX: "auto", 
                                    minHeight:"255px", 
                                    width: "100%", 
                                    pb: 3, 
                                }}
                            >
                                <DragDropContext onDragEnd={onDragEnd} style={{ height: "100%" }}>
                                    <Droppable droppableId="droppable" direction="horizontal">
                                        {(provided) => (
                                        <div ref={provided.innerRef} {...provided.droppableProps} style={{ display: "flex", width: "max-content", height: "100%" }}>
                                            {columns.map((column, index) => (
                                            <Draggable key={index} draggableId={`column-${index}`} index={index}>
                                                {(provided) => (
                                                <div ref={provided.innerRef} {...provided.draggableProps} {...provided.dragHandleProps}>
                                                    <Paper 
                                                        elevation={3} 
                                                        sx={{ backgroundColor: "#f0f0f0" }}
                                                        style={{ 
                                                            minWidth: 200, 
                                                            maxWidth: 200,
                                                            margin: "0 8px",
                                                            marginRight: index === columns.length - 1 ? "8px" : "4px" 
                                                        }}
                                                    >   
                                                        <Box
                                                            sx={{
                                                                backgroundColor: "#555",
                                                                color: "#fff",
                                                                padding: "4px 8px",
                                                                display: "flex",
                                                                justifyContent: "space-between",
                                                                alignItems: "center",
                                                                marginBottom: 1,
                                                                borderTopLeftRadius: "4px",
                                                                borderTopRightRadius: "4px",
                                                                borderBottom: "1px solid #ccc",
                                                            }}
                                                        >
                                                            <Typography variant="h6" component="div">
                                                                Motif {index + 1}
                                                            </Typography>
                                                            <Box sx={{ display: "flex", gap: "3px" }}>
                                                                <Tooltip title="Delete motif.">
                                                                    <IconButton 
                                                                        onClick={() => handleDeleteColumn(index)}
                                                                        style={{ color: "#990000", backgroundColor: "#ff5f57", padding: "1px"}}
                                                                    >
                                                                        <CloseIcon style={{fontSize: "12px"}} />
                                                                    </IconButton>
                                                                </Tooltip>
                                                            </Box>
                                                        </Box>
                                                        <Box style={{ width: "100%" }}>
                                                            {column.map((motif, columnItemIndex) => (
                                                                <Box sx={{backgroundColor: "#ccc", marginBottom: "10px", minHeight: "130px"}}>
                                                                    <Box 
                                                                        sx={{
                                                                            backgroundColor: "#ccc",
                                                                            display: "flex",
                                                                            justifyContent: "space-between",
                                                                            alignItems: "center",
                                                                            marginBottom: 1,
                                                                            padding: "4px 8px",
                                                                            paddingTop: "0px",
                                                                            width: "100%",
                                                                        }}
                                                                    >
                                                                        <FormControlLabel
                                                                            control={<Switch
                                                                                checked={motif.motifType === "polyketide"}
                                                                                onChange={(event) => {
                                                                                    const newColumns = columns.map((col, colIndex) => 
                                                                                        colIndex === index 
                                                                                            ? col.map((item, itemIndex) => 
                                                                                                itemIndex === columnItemIndex 
                                                                                                    ? { ...item, motifType: event.target.checked ? "polyketide" : "peptide" } 
                                                                                                    : item
                                                                                            )
                                                                                            : col
                                                                                    );
                                                                                    setColumns(newColumns);
                                                                                }}
                                                                                color="primary"
                                                                                name={`motif-${index}-${columnItemIndex}`}
                                                                            />}
                                                                            label={motif.motifType === "peptide" ? "Other" : motif.motifType.charAt(0).toUpperCase() + motif.motifType.slice(1)}
                                                                        />
                                                                        <Box sx={{ display: "flex", gap: "3px" }}>
                                                                            {(motif.motifType === "peptide") && (
                                                                                <Tooltip title="Search for ID on PubChem.">
                                                                                    <IconButton 
                                                                                        onClick={() => {
                                                                                            // check if the PubChem ID is valid by seeing if you can parse it as integer
                                                                                            if (isNaN(parseInt(motif.peptideCid))) {
                                                                                                toast.warn("Invalid PubChem ID.");
                                                                                                return;
                                                                                            };
                                                                                            window.open(`https://pubchem.ncbi.nlm.nih.gov/compound/${motif.peptideCid}`, "_blank");
                                                                                        }}
                                                                                        style={{ color: "#985600", backgroundColor: "#febc2e", padding: "1px"}}
                                                                                    >
                                                                                        <SearchIcon style={{fontSize: "12px" }} />
                                                                                    </IconButton>
                                                                                </Tooltip>
                                                                            )}
                                                                            {(column.length > 1) && (
                                                                                <Tooltip title="Delete option.">
                                                                                    <IconButton 
                                                                                        onClick={() => {
                                                                                            const newColumns = columns.map((col, colIndex) => 
                                                                                                colIndex === index 
                                                                                                    ? col.filter((_, itemIndex) => itemIndex !== columnItemIndex)
                                                                                                    : col
                                                                                            );
                                                                                            setColumns(newColumns);
                                                                                        }}
                                                                                        style={{ color: "#990000", backgroundColor: "#ff5f57", padding: "1px"}}
                                                                                    >
                                                                                        <CloseIcon style={{fontSize: "12px"}} />
                                                                                    </IconButton>
                                                                                </Tooltip>
                                                                            )}
                                                                        </Box>
                                                                    </Box>
                                                                    {motif.motifType === "polyketide" && (
                                                                        <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center", justifyContent: "space-between", paddingBottom: "10px" }}>
                                                                            <Box sx={{ width: "100%", paddingLeft: "10px", paddingRight: "2.5px" }}>
                                                                                <FormControl
                                                                                    variant="filled"
                                                                                    style={{ 
                                                                                        width: "100%",
                                                                                        borderRadius: "0px",
                                                                                     }}
                                                                                >
                                                                                    <InputLabel id="demo-simple-select-label">Type</InputLabel>
                                                                                    <Select
                                                                                        labelId="demo-simple-select-label"
                                                                                        id="demo-simple-select"
                                                                                        value={motif.polyketideType}
                                                                                        onChange={(event) => {
                                                                                            const newColumns = columns.map((col, colIndex) => 
                                                                                                colIndex === index 
                                                                                                    ? col.map((item, itemIndex) => 
                                                                                                        itemIndex === columnItemIndex 
                                                                                                            ? { ...item, polyketideType: event.target.value } 
                                                                                                            : item
                                                                                                    )
                                                                                                    : col
                                                                                            );
                                                                                            setColumns(newColumns);
                                                                                        }}
                                                                                        style={{ width: "100%", padding: "4px" }}
                                                                                        MenuProps={{
                                                                                            PaperProps: {
                                                                                                style: {
                                                                                                    maxHeight: 300,
                                                                                                    overflowY: "auto",
                                                                                                },
                                                                                            },   
                                                                                        }}
                                                                                    >
                                                                                        <MenuItem value="Any">Any</MenuItem>
                                                                                        <MenuItem value="A">A</MenuItem>
                                                                                        <MenuItem value="B">B</MenuItem>
                                                                                        <MenuItem value="C">C</MenuItem>
                                                                                        <MenuItem value="D">D</MenuItem>
                                                                                    </Select>
                                                                                </FormControl>
                                                                            </Box>
                                                                            <Box sx={{ width: "100%", paddingRight: "10px", paddingLeft: "2.5px" }}>
                                                                                <FormControl
                                                                                    variant="filled"
                                                                                    style={{
                                                                                        width: "100%",
                                                                                        borderRadius: "0px",
                                                                                    }}
                                                                                >
                                                                                    <InputLabel id="demo-simple-select-label">Decor</InputLabel>
                                                                                    <Select
                                                                                        labelId="demo-simple-select-label"
                                                                                        id="demo-simple-select"
                                                                                        value={motif.polyketideDecor}
                                                                                        onChange={(event) => {
                                                                                            const newColumns = columns.map((col, colIndex) =>   
                                                                                                colIndex === index
                                                                                                    ? col.map((item, itemIndex) =>
                                                                                                        itemIndex === columnItemIndex
                                                                                                            ? { ...item, polyketideDecor: event.target.value }
                                                                                                            : item
                                                                                                    )
                                                                                                    : col
                                                                                            );
                                                                                            setColumns(newColumns);
                                                                                        }}
                                                                                        style={{ width: "100%", padding: "4px" }}
                                                                                        MenuProps={{
                                                                                            PaperProps: {
                                                                                                style: {
                                                                                                    maxHeight: 300,
                                                                                                    overflowY: "auto",
                                                                                                },
                                                                                            },   
                                                                                        }}
                                                                                    >
                                                                                        <MenuItem value="Any">Any</MenuItem>
                                                                                        <MenuItem value="1">1</MenuItem>
                                                                                        <MenuItem value="2">2</MenuItem>
                                                                                        <MenuItem value="3">3</MenuItem>
                                                                                        <MenuItem value="4">4</MenuItem>
                                                                                        <MenuItem value="5">5</MenuItem>
                                                                                        <MenuItem value="6">6</MenuItem>
                                                                                        <MenuItem value="7">7</MenuItem>
                                                                                        <MenuItem value="8">8</MenuItem>
                                                                                        <MenuItem value="9">9</MenuItem>
                                                                                        <MenuItem value="10">10</MenuItem>
                                                                                        <MenuItem value="11">11</MenuItem>
                                                                                        <MenuItem value="12">12</MenuItem>
                                                                                        <MenuItem value="13">13</MenuItem>
                                                                                        <MenuItem value="14">14</MenuItem>
                                                                                        <MenuItem value="15">15</MenuItem>
                                                                                    </Select>
                                                                                </FormControl>
                                                                            </Box>
                                                                        </Box>
                                                                    )}
                                                                    {motif.motifType === "peptide" && (
                                                                        <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center", justifyContent: "space-between", paddingBottom: "6px" }}>
                                                                            <Box sx={{ width: "100%", paddingRight: "10px", paddingLeft: "10px" }}>
                                                                                <TextField
                                                                                    value={motif.peptideCid} 
                                                                                    defaultValue={"Any"}
                                                                                    onChange={(event) => {
                                                                                        const newColumns = columns.map((col, colIndex) =>
                                                                                            colIndex === index
                                                                                                ? col.map((item, itemIndex) =>
                                                                                                    itemIndex === columnItemIndex
                                                                                                        ? { ...item, peptideCid: event.target.value }
                                                                                                        : item
                                                                                                )
                                                                                                : col
                                                                                        );
                                                                                        setColumns(newColumns);
                                                                                    }}
                                                                                    label="PubChem ID"
                                                                                    style={{ width: "100%", padding: "4px" }}
                                                                                />
                                                                            </Box>
                                                                        </Box>
                                                                    )}
                                                                </Box>
                                                            ))}
                                                        </Box>
                                                        <Tooltip title="Add option">
                                                            <Box
                                                                type="button"
                                                                sx={{
                                                                    backgroundColor: "#1976d2",
                                                                    color: "white",
                                                                    height: "35px",
                                                                    display: "flex",
                                                                    alignItems: "center",
                                                                    justifyContent: "center",
                                                                    marginBottom: 1,
                                                                    borderBottomLeftRadius: "4px",
                                                                    borderBottomRightRadius: "4px",
                                                                    cursor: "pointer",
                                                                    width: "100%",
                                                                }}
                                                                onClick={addMotifToColumn(index)}
                                                            >
                                                                <AddIcon style={{ alignSelf: "center" }} />
                                                            </Box>
                                                        </Tooltip>
                                                    </Paper>
                                                </div>
                                                )}
                                            </Draggable>
                                            ))}
                                            {provided.placeholder}
                                        </div>
                                        )}
                                    </Droppable>
                                </DragDropContext>
                            </Box>
                        ) : (
                            <Typography variant="body1" sx={{ mb: 2, paddingLeft: "10px" }}>
                                {`Query is empty.`}
                            </Typography>
                        )}
                    </Box>
                </Box>

                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginTop: "20px" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="Choose between pairwise matching or exact pattern matching.">
                            <Typography variant="h6" gutterBottom>
                                Query type
                            </Typography>
                        </TypographyTooltip>
                    </Box>
                </Box>
                <Divider sx={{ mb: 2 }} />
                <Box sx={{ display: "flex", gap: 2}}>
                    <Box sx={{ 
                        flexGrow: 1, 
                        flex: "1 1 50%",
                        padding: "10px"
                    }}>
                        <RadioGroup
                            value={queryType}
                            onChange={(event) => setQueryType(event.target.value)}
                            sx={{ mb: 2 }}
                        >
                            <FormControlLabel
                                value="match"
                                control={<Radio />}
                                label="Pairwise matching"
                            />
                            <FormControlLabel
                                value="query"
                                control={<Radio />}
                                label="Pattern matching"
                            />
                        </RadioGroup>
                    </Box>
                    <Box sx={{ 
                        flexGrow: 1, 
                        flex: "1 1 50%",
                        padding: "10px"
                    }}>
                        {queryType === "match" && (
                            <RadioGroup
                                value={alignmentType}
                                onChange={(event) => setAlignmentType(event.target.value)}
                                sx={{ mb: 2 }}
                            >
                                <FormControlLabel
                                    value="global"
                                    control={<Radio />}
                                    label="Global alignment strategy"
                                />
                                <FormControlLabel
                                    value="local"
                                    control={<Radio />}
                                    label="Local alignment strategy"
                                />
                            </RadioGroup>
                        )}
                        {queryType === "match" && (
                            <Box>
                                <TextField
                                    value={gapPenalty}
                                    onChange={(event) => setGapPenalty(event.target.value)}
                                    label="Gap penalty"
                                    margin="normal"
                                    variant="outlined"
                                    type="number"
                                    sx={{ width: "100%", mb: 2 }}
                                />
                                <TextField
                                    value={endGapPenalty}
                                    onChange={(event) => setEndGapPenalty(event.target.value)}
                                    label="End gap penalty"
                                    margin="normal"
                                    variant="outlined"
                                    type="number"
                                    sx={{ width: "100%", mb: 2 }}
                                />
                            </Box>
                        )}
                        {queryType === "query" && (
                            <Box>
                                <FormControlLabel
                                    control={<Switch 
                                        checked={queryHasLeadingModules}
                                        onChange={(event) => setQueryHasLeadingModules(event.target.checked)}
                                        color="primary"
                                        name="queryHasLeadingModules"
                                    />}
                                    label="Query has leading modules"
                                />
                                <FormControlLabel
                                    control={<Switch 
                                        checked={queryHasTrailingModules}
                                        onChange={(event) => setQueryHasTrailingModules(event.target.checked)}
                                        color="primary"
                                        name="queryHasTrailingModules"
                                    />}
                                    label="Query has trailing modules"
                                />
                            </Box>
                        )}
                    </Box>
                </Box>

                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginTop: "20px" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="Select filters to filter results on. If no filters are selected, all results will be shown.">
                            <Typography variant="h6" gutterBottom>
                                Filters
                            </Typography>
                        </TypographyTooltip>
                    </Box>
                </Box>
                <Divider sx={{ mb: 2 }} />
                <Box sx={{ display: "flex", gap: 2 }}>
                    <MultiSelect 
                        title="Select bioactivity labels"
                        labels={allBioactivityLabels}
                        selectedLabels={selectedBioactivityLabels}
                        setSelectedLabels={setSelectedBioactivityLabels}
                    />
                    <MultiSelect 
                        title="Select organisms"
                        labels={allOrganismLabels}
                        selectedLabels={selectedOrganismLabels}
                        setSelectedLabels={setSelectedOrganismLabels}
                    />
                </Box>

                <Box sx={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginTop: "20px" }}>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TypographyTooltip title="Further define the scope of the query.">
                            <Typography variant="h6" gutterBottom>
                                Query space
                            </Typography>
                        </TypographyTooltip>
                    </Box>
                </Box>
                <Divider sx={{ mb: 2 }} />
                <Box sx={{ display: "flex", gap: 2 }}>
                    <Box sx={{ 
                        flexGrow: 1, 
                        flex: "1 1 50%",
                        padding: "10px"
                    }}>
                        <Box>
                            <FormControlLabel
                                control={<Switch 
                                    checked={queryAgainstMolecules}
                                    onChange={(event) => setQueryAgainstMolecules(event.target.checked)}
                                    color="primary"
                                    name="queryAgainstMolecules"
                                />}
                                label="Query against molecules"
                            />
                            <FormControlLabel
                                control={<Switch 
                                    checked={queryAgainstProtoclusters}
                                    onChange={(event) => setQueryAgainstProtoclusters(event.target.checked)}
                                    color="primary"
                                    name="queryAgainstProtoclusters"
                                />}
                                label="Query against protoclusters"
                            />
                        </Box>
                    </Box>
                    <Box sx={{ 
                        flexGrow: 1, 
                        flex: "1 1 50%",
                        padding: "10px"
                    }}>
                        <TextField
                            value={maxNumMatches}
                            onChange={(event) => setMaxNumMatches(event.target.value)}
                            label="Max shown matches"
                            margin="normal"
                            variant="outlined"
                            type="number"
                            sx={{ width: "100%", mb: 2 }}
                        />
                        <TextField
                            value={minMatchLength}
                            onChange={(event) => setMinMatchLength(event.target.value)}
                            label="Min length motif code"
                            margin="normal"
                            variant="outlined"
                            type="number"
                            sx={{ width: "100%", mb: 2 }}
                        />
                        <TextField
                            value={maxMatchLength}
                            onChange={(event) => setMaxMatchLength(event.target.value)}
                            label="Max length motif code"
                            margin="normal"
                            variant="outlined"
                            type="number"
                            sx={{ width: "100%", mb: 2 }}
                        />
                    </Box>
                </Box>
                
                <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center", marginTop: "20px" }}>
                    <Button
                        variant="contained"
                        color="primary"
                        onClick={() => submit({
                            query: columns,
                            queryType: queryType,
                            alignmentType: alignmentType,
                            gapPenalty: gapPenalty,
                            endGapPenalty: endGapPenalty,
                            queryHasLeadingModules: queryHasLeadingModules,
                            queryHasTrailingModules: queryHasTrailingModules,
                            selectedBioactivityLabels: selectedBioactivityLabels,
                            selectedOrganismLabels: selectedOrganismLabels,
                            maxNumMatches: maxNumMatches,
                            queryAgainstMolecules: queryAgainstMolecules,
                            queryAgainstProtoclusters: queryAgainstProtoclusters,
                            minMatchLength: minMatchLength,
                            maxMatchLength: maxMatchLength,
                        })}
                        style={{ width: "150px" }}
                    >
                        Submit query
                    </Button>
                </Box>

            </Box>
        </Paper>
    );
};

export default QueryForm;