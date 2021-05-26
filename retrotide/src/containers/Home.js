import React from 'react';
import {
  Link
} from "react-router-dom";
import Button from '../components/Button';

// all of this needs to move to presentational later but we'll see
// if we end up using that paradigm
function Home() {
	return (
    <div className='home form'>
      <h3>Welcome</h3>
      <div>
        <p><b>Retrotide</b> is a retrobiosynthesis tool for PKS design by researchers and institutions. For more information see the About page or
        select an option below or in the navigation tabs: </p>
        <table className='welcomeTable'>
          <tbody>
            <tr>
              <td>Design PKS</td>
              <td>Start here to begin a design process for PKS modules that will generate a target molecule</td>
              <td><Link to='/Design'><Button>Click Here >></Button></Link></td>
            </tr>
            <tr>
              <td>Architecture Search</td>
              <td>(in progress)</td>
              <td><Link to='/Search'><Button>Click Here >></Button></Link></td>
            </tr>
            <tr>
              <td>Saved Designs</td>
              <td>Existing users can login to view previously-run design projects</td>
              <td><Link to='/SavedResults'><Button>Click Here >></Button></Link></td>
            </tr>
            <tr>
              <td>About</td>
              <td>More information about Retrotide including citation information</td>
              <td><Link to='/About'><Button>Click Here >></Button></Link></td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>
  )
};

export default Home;