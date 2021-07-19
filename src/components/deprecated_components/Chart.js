import React from 'react';
import { LineChart, Line, XAxis, YAxis, Label, ResponsiveContainer } from 'recharts';
import { makeStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import Info from '@material-ui/icons/Info';
import IconButton from '@material-ui/core/IconButton';
import Chip from '@material-ui/core/Chip';
import Typography from '@material-ui/core/Typography';

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexWrap: 'wrap',
  },
  chip: {
    margin: theme.spacing(1),
  },
  title: {
    flexGrow:1,
  }
}));

// Generate Sales Data
function createData(time, amount) {
  return { time, amount };
}

const data = [
  createData('00:00', 0),
  createData('03:00', 300),
  createData('06:00', 600),
  createData('09:00', 800),
  createData('12:00', 1500),
  createData('15:00', 2000),
  createData('18:00', 2400),
  createData('21:00', 2400),
  createData('24:00', undefined),
];

export default function Chart() {
  const classes = useStyles();
  return (
    <React.Fragment>
    <div className={classes.root}>
      <Typography component="h2" variant="h6" color="primary" gutterBottom>
        Arbitrary Figure
        <IconButton aria-label="info">
          <Info fontSize="small"/>
        </IconButton>
      </Typography>
      <Chip label="Sample 1" className={classes.chip} />
      <Chip label="Sample 2" className={classes.chip} />
      <Chip label="Sample 3" className={classes.chip} />
      </div>
      <ResponsiveContainer>
        <LineChart
          data={data}
          margin={{
            top: 16,
            right: 16,
            bottom: 0,
            left: 24,
          }}
        >
          <XAxis dataKey="time" />
          <YAxis>
            <Label angle={270} position="left" style={{ textAnchor: 'middle' }}>
              y-Axis
            </Label>
          </YAxis>
          <Line type="monotone" dataKey="amount" stroke="#556CD6" dot={false} />
        </LineChart>
      </ResponsiveContainer>
    </React.Fragment>
  );
}