import React from 'react';
import { Link } from 'react-router-dom';
import { withRouter } from 'react-router';

// navigation link with some logic idenitfying active tab

function NavTab(props) {

  const {children, className, path, location} = props;
  const currentPath = location.pathname;
  const active = (path === currentPath);
  let classFinal = (className ? className + " " : "") + "Navbutton" + (active ? " activeTab" : "")

  return (
    <Link to={ path } className={ classFinal } >
      { children }
    </Link>
  )
}


export default withRouter(NavTab);