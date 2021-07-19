import React from 'react';
import clsx from 'clsx';
import { makeStyles } from '@material-ui/core/styles';
import CssBaseline from '@material-ui/core/CssBaseline';
import Drawer from '@material-ui/core/Drawer';
import AppBar from '@material-ui/core/AppBar';
import Toolbar from '@material-ui/core/Toolbar';
import List from '@material-ui/core/List';
import Typography from '@material-ui/core/Typography';
import Divider from '@material-ui/core/Divider';
import IconButton from '@material-ui/core/IconButton';
import Badge from '@material-ui/core/Badge';
import Container from '@material-ui/core/Container';
import Grid from '@material-ui/core/Grid';
import Paper from '@material-ui/core/Paper';
import Link from '@material-ui/core/Link';
import MenuIcon from '@material-ui/icons/Menu';
import ChevronLeftIcon from '@material-ui/icons/ChevronLeft';
import NotificationsIcon from '@material-ui/icons/Notifications';
import { mainListItems } from './listItems';
import Header from './Header';
import Chart from './Chart';
import Button from '@material-ui/core/Button';
import Info from '@material-ui/icons/Info';
import TextField from '@material-ui/core/TextField';
import Card from '@material-ui/core/Card';
import CardActions from '@material-ui/core/CardActions';
import Select from '@material-ui/core/Select';
import FormControl from '@material-ui/core/FormControl';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import MenuItem from '@material-ui/core/MenuItem';
import InputLabel from '@material-ui/core/InputLabel';
import FigureSettings from './FigureSettings';

function MadeWithLove() {
  return (
    <Typography variant="body2" color="textSecondary" align="center">
      {'Built by Joshua Wu of the '}
      <Link color="inherit" href="http://www.jasonspencelab.com/">
        Jason Spence Lab
      </Link>
      {'.'}
    </Typography>
  );
}

const drawerWidth = 240;

const analysisFields = [
  {
    id:'n_neighbors',
    label:'Nearest Neighbors',
    defaultValue:'15',
    helperText:'Size of local neighborhood'
  },
  {
    id:'n_pcs',
    label:'Principle Components',
    defaultValue:'11',
    helperText:'# of principle components used in construction of neighborhood graph'
  },
  {
    id:'spread',
    label:'Spread',
    defaultValue:'1',
    helperText:'How clumped embedded points are'
  },
  {
    id:'min_dist',
    label:'Min Distance',
    defaultValue:'0.4',
    helperText:'Minimum distance between points on the cluster graph'
  },
  {
    id:'resolution',
    label:'Resolution',
    defaultValue:'0.1',
    helperText:'Number of clusters identified'
  }
]

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
  },
  toolbar: {
    paddingRight: 24, // keep right padding when drawer closed
  },
  toolbarIcon: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'flex-end',
    padding: '0 8px',
    ...theme.mixins.toolbar,
  },
  appBar: {
    zIndex: theme.zIndex.drawer + 1,
    transition: theme.transitions.create(['width', 'margin'], {
      easing: theme.transitions.easing.sharp,
      duration: theme.transitions.duration.leavingScreen,
    }),
    borderBottom: `0px solid ${theme.palette.divider}`,
    background: 'linear-gradient(45deg, #041E42 30%, #041E42 90%)',
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
  title: {
    flexGrow: 1,
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
  appBarSpacer: theme.mixins.toolbar,
  content: {
    flexGrow: 1,
    height: '100vh',
    overflow: 'auto',
  },
  container: {
    paddingTop: theme.spacing(4),
    paddingBottom: theme.spacing(4),
  },
  paper: {
    padding: theme.spacing(2),
    display: 'flex',
    overflow: 'auto',
    flexDirection: 'column',
  },
  fixedHeight: {
    height: 500,
  },
  link: {
    margin: theme.spacing(1, 1.5),
  },
  submit: {
        margin: theme.spacing(0,1,1),
        background: 'linear-gradient(45deg, #041E42 30%, #041E42 90%)',
  },
}));

export default function Dashboard() {
  const classes = useStyles();
  const [open, setOpen] = React.useState(true);
  const handleDrawerOpen = () => {
    setOpen(true);
  };
  const handleDrawerClose = () => {
    setOpen(false);
  };
  const fixedHeightPaper = clsx(classes.paper, classes.fixedHeight);

  return (
    <div className={classes.root}>
      <Header inputPage={false}/>

      <main className={classes.content}>
        <div className={classes.appBarSpacer} />
        <Container maxWidth="lg" className={classes.container}>
          <Grid container spacing={3}>
            {/* Chart */}
            <Grid item xs={12} md={8} lg={9}>
              <Paper className={fixedHeightPaper} display="flex">
                <Chart />
              </Paper>
            </Grid>
            {/* Recent Deposits */}
            <Grid item xs={12} md={4} lg={3}>
              <Paper className={fixedHeightPaper}>
                <FigureSettings/>
              </Paper>
            </Grid>
            {/* Recent Orders */}
            <Grid item xs={12} md={8} lg={9}>
              <Paper className={classes.paper}>
                <Typography component="h1" variant="h6" color="inherit" noWrap className={classes.title}>
                  Analysis Parameters
                  <IconButton className={classes.button} aria-label="info">
                    <Info fontSize="small"/>
                  </IconButton>
                </Typography>
                <Grid container spacing={1}>
                      {analysisFields.map(field =>(
                        <Grid item xs={6} sm={4}>
                        <TextField
                            id={field.id}
                            label={field.label}
                            style={{ margin: 8 }}
                            defaultValue={field.defaultValue}
                            //helperText={field.helperText}
                            fullWidth
                            margin="normal"
                            InputLabelProps={{
                              shrink: true,
                            }}
                        />
                        </Grid>
                    ))}
                  </Grid>
                
              </Paper>
            </Grid>
            <Grid container xs={12} md={4} lg={3} direction="column" alignItems="center" justify="center">
              <Grid item>
                <Button
                      type="submit"
                      variant="contained"
                      color="secondary"
                      size="large"
                      className={classes.submit}
                      >
                    Update Figure
                </Button>
              </Grid>
            </Grid>
          </Grid>
        </Container>
        <MadeWithLove />
      </main>
    </div>
  );
}