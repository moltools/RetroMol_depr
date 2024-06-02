import React, { useEffect, useState } from "react";
import { Box, Button, Divider, Grid, IconButton, Paper, Tooltip, Typography } from "@mui/material";
import { DragDropContext, Droppable, Draggable } from "react-beautiful-dnd";
import TypographyTooltip from "../common/TypographyTooltip";
import RefreshIcon from "@mui/icons-material/Refresh";
import DeleteIcon from "@mui/icons-material/Delete";
import CloseIcon from "@mui/icons-material/Close";
import AddIcon from "@mui/icons-material/Add";

const QueryForm = ({ results, selectedResultIndex }) => {
    const [columns, setColumns] = useState([]);

    const handleRefresh = () => {
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
        const newColumns = [...columns, []];
        setColumns(newColumns);
    };

    const addMotifToColumn = (index) => () => {
        const newColumns = Array.from(columns);
        newColumns[index] = [...newColumns[index], "New"];
        setColumns(newColumns);
    };

    useEffect(() => {
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
                            <Box sx={{ overflowX: "auto", minHeight:"255px", width: "100%", pb: 3, }}>
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
                                                            // padding: 16, 
                                                            minWidth: 150, 
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
                                                            }}
                                                        >
                                                            <Typography variant="h6" component="div">
                                                                Motif {index + 1}
                                                            </Typography>
                                                            <Box sx={{ display: "flex", gap: "3px" }}>
                                                                <Tooltip title="Delete motif">
                                                                    <IconButton 
                                                                        onClick={() => handleDeleteColumn(index)}
                                                                        style={{ color: "#990000", backgroundColor: "#ff5f57", padding: "1px"}}
                                                                    >
                                                                        <CloseIcon style={{fontSize: "12px"}} />
                                                                    </IconButton>
                                                                </Tooltip>
                                                            </Box>
                                                        </Box>
                                                        <Box style={{ padding: "16px" }}>
                                                            {column.map((item, idx) => (
                                                                <Typography key={idx} variant="body1">
                                                                    {item}
                                                                </Typography>
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
            </Box>
        </Paper>
    );
};

export default QueryForm;