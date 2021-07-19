import React from 'react';
import Grid from '@material-ui/core/Grid';
import { withStyles } from '@material-ui/core/styles';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import FormControl from '@material-ui/core/FormControl';
import IconButton from '@material-ui/core/IconButton';
import Info from '@material-ui/icons/Info';
import InputLabel from '@material-ui/core/InputLabel';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField'

const styles = theme => ({
  root: {
    display: 'flex',
  },
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120,
  },
  title: {
    flexGrow: 1,
  },
})

class FigureSettings extends React.Component {
  constructor(props) {
    super(props);
  }
  render() {
    const {classes} = this.props;
    return(
      <div>
        <Typography component="h1" variant="h6" color="inherit" noWrap className={classes.title}>
          Figure Settings
          <IconButton className={classes.button} aria-label="info">
            <Info fontSize="small"/>
          </IconButton>
        </Typography>
        <Grid container spacing={1}>
          <Grid item xs={6} sm={6} md={6} lg={6}>
          <FormControl className={classes.formControl}>
            <InputLabel htmlFor="age-simple">Color Palette</InputLabel>
            <Select
              // onChange={handleChange}
              inputProps={{
                name: 'age',
                id: 'age-simple',
              }}
            >
              <MenuItem value={10}>Ten</MenuItem>
              <MenuItem value={20}>Twenty</MenuItem>
              <MenuItem value={30}>Thirty</MenuItem>
            </Select>
          </FormControl>
          <FormControl className={classes.formControl}>
            <InputLabel htmlFor="age-simple">Color Map</InputLabel>
            <Select
              // onChange={handleChange}
              inputProps={{
                name: 'age',
                id: 'age-simple',
              }}
            >
              <MenuItem value={10}>Ten</MenuItem>
              <MenuItem value={20}>Twenty</MenuItem>
              <MenuItem value={30}>Thirty</MenuItem>
            </Select>
          </FormControl>

          <TextField
              id={1}
              label="Dot Size"
              style={{ margin: 8 }}
              defaultValue={5}
              //helperText={field.helperText}
              fullWidth
              margin="normal"
              InputLabelProps={{
                shrink: true,
              }}
          />
          </Grid>
        </Grid>
      </div>
    )
  }
}

export default withStyles(styles)(FigureSettings);