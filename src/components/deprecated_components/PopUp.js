import React from 'react';
import PropTypes from 'prop-types';
import { withStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import Avatar from '@material-ui/core/Avatar';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemAvatar from '@material-ui/core/ListItemAvatar';
import ListItemText from '@material-ui/core/ListItemText';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogContent from '@material-ui/core/DialogContent';
import Dialog from '@material-ui/core/Dialog';
import PersonIcon from '@material-ui/icons/Person';
import AddIcon from '@material-ui/icons/Add';
import Typography from '@material-ui/core/Typography';
import { blue } from '@material-ui/core/colors';

const emails = ['username@gmail.com', 'user02@gmail.com'];
const styles = theme => ({
  avatar: {
    backgroundColor: blue[100],
    color: blue[600],
  },
});


class SimpleDialogDemo extends React.Component {
  constructor(props) {
    super(props);
    this.handleClickOpen = this.handleClickOpen.bind(this);
    this.handleClose = this.handleClose.bind(this);
    this.state = {
      open:false,
    }
  }

  handleClickOpen(e) {
    this.setState({open:true});
  }

  handleClose(e) {
    this.setState({open:false});
  };
  render() {
  return (
    <div>
      <Button variant="outlined" color="primary" onClick={this.handleClickOpen}>
        Open simple dialog
      </Button>
      <Dialog onClose={this.handleClose} aria-labelledby="simple-dialog-title" open={this.state.open}>
        <DialogTitle id="simple-dialog-title">Set backup account</DialogTitle>
        <DialogContent>

          sdlfkj
    
        </DialogContent>
        <DialogTitle id="simple-dialog-title">Set backup account</DialogTitle>
        <DialogContent>

          sdlfkj
    
        </DialogContent>
      </Dialog>
    </div>
  )};
}

export default withStyles(styles)(SimpleDialogDemo)