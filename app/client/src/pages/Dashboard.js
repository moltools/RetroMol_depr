import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";
import { Link } from "react-router-dom";
import { 
    AppBar, 
    Box, 
    Container,
    Divider,
    Drawer, 
    Grid,
    IconButton, 
    List, 
    ListItem, 
    ListItemIcon, 
    ListItemText, 
    Toolbar, 
    Typography 
} from "@mui/material";
import { 
    BugReport as BugReportIcon,
    Home as HomeIcon, 
    Info as InfoIcon, 
    Menu as MenuIcon
} from "@mui/icons-material";

import InputForm from "../components/dashboard/InputForm";
import LoadingOverlay from "../components/common/LoadingOverlay";
import Page from "../components/common/Page";
import Status from "../components/common/Status";

function Dashboard () {
    const [isBusy, setIsBusy] = useState(false);
    const [isDrawerOpen, setIsDrawerOpen] = useState(false);
    const [version, setVersion] = useState("0.0.0");
    const [statusServer, setStatusServer] = useState(true);
    const [statusDatabase, setStatusDatabase] = useState(false);
    const [parseResults, setParseResults] = useState([]);

    const toggleDrawer = () => {
        setIsDrawerOpen(!isDrawerOpen);
    };

    const fetchVersion = async () => {
        try {
            const response = await fetch("/api/fetch_version");
            const data = await response.json();
            if (data.status === "success") {
                setVersion(data.payload.version);
            };
        } catch (error) {
            console.error(error);
        };
    };

    const checkStatusServer = async () => {
        try {
            const response = await fetch("/api/fetch_server_status");
            const data = await response.json();
            if (data.status === "success") {
                setStatusServer(true);
            };
        } catch (error) {
            console.error(error);
        };
    };

    const checkStatusDatabase = async () => {
        try {
            const response = await fetch("/api/fetch_database_status");
            const data = await response.json();
            if (data.status === "success") {
                setStatusDatabase(true);
            };
        } catch (error) {
            console.error(error);
        };
    };

    useEffect(() => {
        fetchVersion();

        if (isDrawerOpen) {
            checkStatusServer();
            checkStatusDatabase();
        }
    }, [isDrawerOpen]);

    const submit = async (data) => {
        setIsBusy(true);

        try {
            const response = await fetch("/api/parse_submission", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify(data)
            });

            const result = await response.json();
            if (result.status === "success") {
                toast.success(result.message);
            } else if (result.status === "warning") {
                toast.warn(result.message);
            } else {
                toast.error(result.message);
            };
            
        } catch (error) {
            console.error(error);
            toast.error("Something went wrong. Please file a bug report.");
        };

        setIsBusy(false);
    };

    return (
        <Page sx={{ display: "flex", width: "100%", flexDirection: "column" }}>
            <AppBar position="fixed">
                <Toolbar>
                    <IconButton
                        color="inherit"
                        aria-label="open drawer"
                        edge="start"
                        onClick={toggleDrawer}
                        sx={{ mr: 2 }}
                    >
                        <MenuIcon />
                    </IconButton>
                    <Typography variant="h6" noWrap>
                        {`RetroMol (${version})`}
                    </Typography>
                </Toolbar>
            </AppBar>
            <Drawer
                variant="temporary"
                open={isDrawerOpen}
                onClose={toggleDrawer}
                ModalProps={{ keepMounted: true }}
                sx={{ '& .MuiDrawer-paper': { width: 240 } }}
            >
                <List>
                    <ListItem 
                        type="button" 
                        component={Link} 
                        to="/"
                        sx={{ textDecoration: "none", color: "#222" }}
                    >
                        <ListItemIcon>
                            <HomeIcon sx={{ color: "#222" }} />
                        </ListItemIcon>
                        <ListItemText primary="Home" />
                    </ListItem>
                    <ListItem 
                        type="button" 
                        component={Link} 
                        to="/about"
                        sx={{ textDecoration: "none", color: "#222" }}
                    >
                        <ListItemIcon>
                            <InfoIcon sx={{ color: "#222" }} />
                        </ListItemIcon>
                        <ListItemText primary="About" />
                    </ListItem>
                    <ListItem 
                        type="button" 
                        component="a" 
                        href="https://github.com/moltools/RetroMol/issues" 
                        target="_blank"
                        sx={{ textDecoration: "none", color: "#222" }}
                    >
                        <ListItemIcon>
                            <BugReportIcon sx={{ color: "#222" }} />
                        </ListItemIcon>
                        <ListItemText primary="Report Bug" />
                    </ListItem>
                    <Divider sx={{ mt: 2, mb: 2 }} />
                    <Status statusName="Server" status={statusServer} />
                    <Status statusName="Database" status={statusDatabase} />
                </List>
            </Drawer>
            <Box
                component="main"
                sx={{ flexGrow: 1, p: 3 }}
            >
                <Toolbar />
                <Container>
                    {isBusy && <LoadingOverlay />}
                    <Grid container spacing={3}>
                        <Grid item xs={12}>
                            <InputForm onSubmit={(data) => submit(data)} />
                        </Grid>
                    </Grid>
                </Container>
            </Box>
        </Page>
    );
};

export default Dashboard;