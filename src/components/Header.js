import React from 'react';
import clsx from 'clsx';
import AppBar from '@material-ui/core/AppBar';
import Button from '@material-ui/core/Button';
import Toolbar from '@material-ui/core/Toolbar';
import Typography from '@material-ui/core/Typography';
import { withStyles } from '@material-ui/core/styles';
import Paper from '@material-ui/core/Paper';
import CssBaseline from '@material-ui/core/CssBaseline';
import IconButton from '@material-ui/core/IconButton';
import MenuIcon from '@material-ui/icons/Menu';
import Drawer from '@material-ui/core/Drawer';
import ChevronLeftIcon from '@material-ui/icons/ChevronLeft';
import Divider from '@material-ui/core/Divider';
import List from '@material-ui/core/List';
// import { mainListItems } from './listItems';
import logo from './u-m_logo-hex.png'

const drawerWidth = 240;

const styles = theme => ({
    root: {
      display: 'flex',
    },
  	appBar: {
        zIndex: theme.zIndex.drawer + 1,
        transition: theme.transitions.create(['width', 'margin'], {
            easing: theme.transitions.easing.sharp,
            duration: theme.transitions.duration.leavingScreen,
        }),
        	borderBottom: `0px solid ${theme.palette.divider}`,
        	background: 'linear-gradient(45deg, #00274C 30%, #00274C 90%)',
        color: 'white'
  	},
    appBarShift: {
        marginLeft: drawerWidth,
        width: `calc(100% - ${drawerWidth}px)`,
        transition: theme.transitions.create(['width', 'margin'], {
            easing: theme.transitions.easing.sharp,
            duration: theme.transitions.duration.enteringScreen,
        }),
    },
    menuButton: {
        marginRight: 36,
    },
    menuButtonHidden: {
        display: 'none',
    },
	toolbar: {
        paddingRight:24,
    	// flexWrap: 'wrap',
  	},
  	toolbarTitle: {
    	flexGrow: 1,
        color: "#FFCB05",
    },
    toolbarIcon: {
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'flex-end',
        padding: '0 8px',
        ...theme.mixins.toolbar,
    },
	link: {
    	margin: theme.spacing(1, 1.5),
        color:"#FFCB05"
    },
    drawerPaper: {
        position: 'relative',
        whiteSpace: 'nowrap',
        width: drawerWidth,
        transition: theme.transitions.create('width', {
          easing: theme.transitions.easing.sharp,
          duration: theme.transitions.duration.enteringScreen,
        }),
    },
    drawerPaperClose: {
        overflowX: 'hidden',
        transition: theme.transitions.create('width', {
          easing: theme.transitions.easing.sharp,
          duration: theme.transitions.duration.leavingScreen,
        }),
        width: theme.spacing(7),
        [theme.breakpoints.up('sm')]: {
          width: theme.spacing(9),
        },
    },
    logo: {
        width: '50px',
        height: '50px',
        marginLeft: '-10px',
        marginRight: '10px',
    },
});

class Header extends React.Component {
    constructor(props) {
      super(props);
      this.handleDrawerOpen = this.handleDrawerOpen.bind(this);
      this.handleDrawerClose = this.handleDrawerClose.bind(this);
      this.state = {open:false}
    }

    handleDrawerOpen(e) {
        this.setState({open:true});
    }

    handleDrawerClose(e) {
        this.setState({open:false})
    }
	
    render() {
        const {classes} = this.props;
        const open = this.state.open;
        const inputPage = this.props.inputPage
        var position = "absolute"
        if (inputPage) {
            position = "static"
        }
    	return (
            <div className={classes.root}>
                
          		<AppBar position={position} color="default" elevation={2} className={clsx(classes.appBar, !inputPage && open && classes.appBarShift)} square={true}>
            		<Toolbar className={classes.toolbar}>
                        {!inputPage &&
                            <IconButton
                                edge="start"
                                color="inherit"
                                aria-label="Open drawer"
                                onClick={this.handleDrawerOpen}
                                className={clsx(classes.menuButton, open && classes.menuButtonHidden)}
                            >
                                <MenuIcon />
                            </IconButton>
                        }
                        <img src={logo} className={classes.logo} alt="Logo" />                  		
                        <Typography component="h1" variant="h6" color="inherit" noWrap className={classes.toolbarTitle} noWrap>
                    			Single Cell RNA Sequencing Analysis Tool
                  		</Typography>
                  		<Button href="http://www.jasonspencelab.com/" color="inherit" variant="outlined" className={classes.link}>
                    			Learn More
                  		</Button>
            		</Toolbar>
            	</AppBar>

{/*                {!inputPage &&
                    <Drawer
                        variant="permanent"
                        classes={{
                          paper: clsx(classes.drawerPaper, !open && classes.drawerPaperClose),
                        }}
                        open={open}
                    >
                        <div className={classes.toolbarIcon}>
                            <IconButton onClick={this.handleDrawerClose}>
                                <ChevronLeftIcon />
                            </IconButton>
                        </div>
                        <Divider />
                            <List>{mainListItems}</List>
                        <Divider /> 
                    </Drawer>
                }*/}
            </div>
        )
    }
};

export default withStyles(styles)(Header);