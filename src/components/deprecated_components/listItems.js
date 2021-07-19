import React from 'react';
import ListItem from '@material-ui/core/ListItem';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemText from '@material-ui/core/ListItemText';
import ListSubheader from '@material-ui/core/ListSubheader';
import DashboardIcon from '@material-ui/icons/Dashboard';
import ShoppingCartIcon from '@material-ui/icons/ShoppingCart';
import PeopleIcon from '@material-ui/icons/People';
import BarChartIcon from '@material-ui/icons/BarChart';
import LayersIcon from '@material-ui/icons/Layers';
import AssignmentIcon from '@material-ui/icons/Assignment';
import ScatterPlotIcon from '@material-ui/icons/ScatterPlot';
import SubjectIcon from '@material-ui/icons/Subject';
import ShowChartIcon from '@material-ui/icons/ShowChart';
import DragIndicatorIcon from '@material-ui/icons/DragIndicator';
import MusicNoteIcon from '@material-ui/icons/MusicNote';

export const mainListItems = (
  <div>
    <ListItem button>
      <ListItemIcon>
        <ScatterPlotIcon/>
      </ListItemIcon>
      <ListItemText primary="UMAP" />
    </ListItem>
    <ListItem button>
      <ListItemIcon>
        <MusicNoteIcon/>
      </ListItemIcon>
      <ListItemText primary="Violin Plot" />
    </ListItem>
    <ListItem button>
      <ListItemIcon>
        <DragIndicatorIcon/>
      </ListItemIcon>
      <ListItemText primary="Dot Plot" />
    </ListItem>
    <ListItem button>
      <ListItemIcon>
        <ShowChartIcon />
      </ListItemIcon>
      <ListItemText primary="PCA" />
    </ListItem>
    <ListItem button>
      <ListItemIcon>
        <SubjectIcon/>
      </ListItemIcon>
      <ListItemText primary="Summary" />
    </ListItem>
  </div>
);

// export const secondaryListItems = (
//   <div>
//     <ListSubheader inset>Saved reports</ListSubheader>
//     <ListItem button>
//       <ListItemIcon>
//         <AssignmentIcon />
//       </ListItemIcon>
//       <ListItemText primary="Current month" />
//     </ListItem>
//     <ListItem button>
//       <ListItemIcon>
//         <AssignmentIcon />
//       </ListItemIcon>
//       <ListItemText primary="Last quarter" />
//     </ListItem>
//     <ListItem button>
//       <ListItemIcon>
//         <AssignmentIcon />
//       </ListItemIcon>
//       <ListItemText primary="Year-end sale" />
//     </ListItem>
//   </div>
// );