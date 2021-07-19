import React from 'react';
import Grid from '@material-ui/core/Grid';
import { withStyles } from '@material-ui/core/styles';
import Box from '@material-ui/core/Box';
import Button from '@material-ui/core/Button';

import Header from './Header';
import Footer from './Footer';
import BasicForm from './BasicForm';
import FormControlLabel from '@material-ui/core/ExpansionPanelDetails';
import Checkbox from '@material-ui/core/Checkbox';
import background from './background.png'

const styles = theme => ({
    '@global': {
    body: {
      	backgroundColor: theme.palette.common.white,
      	backgroundImage: `url(${background})`,
        backgroundRepeat: 'no-repeat',
        backgroundSize: 'cover',
	    backgroundPosition: 'center',
	    // fontFamily: ['Literata', 'serif'].join(','),
    },},
    containerGrid: {
        margin: 0,
        width: '100%',
    },
    headerGrid: {
        padding: "0px",
        marginBottom:"10px",
    },
    paramGrid: {
        padding: "4px",
        marginBottom:"5px",
    }
});

class FormPage extends React.Component {
    constructor(props) {
        super(props);
        this.handleOnClick = this.handleOnClick.bind(this);
        this.state = {
            isInputPage: true,
        }
    }

    handleOnClick(e) {
        this.setState({isInputPage:false})
    }
    
    render() {
        const {classes} = this.props;
        return (
            <React.Fragment>
            <Grid container wrap={"wrap"} component="main" xs={12} className={classes.containerGrid}>
        	    <Grid item sm={12} md={12} className={classes.headerGrid}>
    	    	    <Header inputPage={this.state.isInputPage}/>
    	        </Grid>
{                <Grid item sm={12} md={12} className={classes.paramGrid}>
    	            <BasicForm/>
                </Grid>}
            </Grid>
            
            <Footer />
            </React.Fragment>
    )};
}

export default withStyles(styles)(FormPage);