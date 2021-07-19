import React from 'react';
import Typography from '@material-ui/core/Typography';
import Link from '@material-ui/core/Link';
import { makeStyles } from '@material-ui/core/styles';

const useStyles = makeStyles(theme => ({
  text: {
  	color:'#FFCB05',
  },
}));

export default function Footer() {
	const classes = useStyles();
  return (
    <Typography variant="body2" color="colorTextSecondary" align="center" className={classes.text}>
      {'Built by Joshua Wu of the '}
      <Link color="inherit" href="http://www.jasonspencelab.com/">
        Jason Spence Lab
      </Link>
      {'.'}
    </Typography>
  );
}